"""
======================
Plasmid Replication
======================

Adapted from chromosome replication.
Performs initiation, elongation, and termination of active plasmid molecules
that replicate independently of the chromosome. (ColE1 - pBR322)

Plasmid replication is initiated asynchronously depending on resource availability
at individual plasmid molecules with no copy number control yet. Second, replication forks
are elongated unidirectionally up to the maximal expected elongation rate, dNTP resource
limitations, and template strand sequence but elongation does not take into account the
action of topoisomerases or the enzymes in the replisome. Finally, replication forks terminate
once they reach the end of their template strand producing fully replicated plasmid molecules
that remain separate from the chromosome and each other.

# TODO: Implement copy number control
"""

import numpy as np

from ecoli.library.schema import (
    numpy_schema,
    counts,
    attrs,
    bulk_name_to_idx,
    listener_schema,
)

from wholecell.utils import units
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease

from ecoli.processes.registries import topology_registry
from ecoli.processes.partition import PartitionedProcess


# Register default topology for this process, associating it with process name
NAME = "ecoli-plasmid-replication"
TOPOLOGY = {
    "bulk": ("bulk",),
    "plasmid_active_replisomes": ("unique", "plasmid_active_replisome"),
    "oriVs": ("unique", "oriV"),
    "plasmid_domains": ("unique", "plasmid_domain"),
    "full_plasmids": ("unique", "full_plasmid"),
    "listeners": ("listeners",),
    "environment": ("environment",),
    "timestep": ("timestep",),
}
topology_registry.register(NAME, TOPOLOGY)


class PlasmidReplication(PartitionedProcess):
    """Plasmid Replication PartitionedProcess"""

    name = NAME
    topology = TOPOLOGY
    defaults = {
        "get_dna_critical_mass": lambda doubling_time: units.Unum,
        "criticalInitiationMass": 975 * units.fg,
        "nutrientToDoublingTime": {},
        "replichore_lengths": np.array([]),
        "sequences": np.array([]),
        "polymerized_dntp_weights": [],
        "replication_coordinate": np.array([]),
        "D_period": np.array([]),
        "replisome_protein_mass": 0,
        "no_child_place_holder": -1,
        "basal_elongation_rate": 967,
        "make_elongation_rates": (
            lambda random, replisomes, base, time_step: units.Unum
        ),
        "mechanistic_replisome": True,
        # molecules
        "replisome_trimers_subunits": [],
        "replisome_monomers_subunits": [],
        "dntps": [],
        "ppi": [],
        # random seed
        "seed": 0,
        "emit_unique": False,
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        # Load parameters
        self.get_dna_critical_mass = self.parameters["get_dna_critical_mass"]
        self.criticalInitiationMass = self.parameters["criticalInitiationMass"]
        self.nutrientToDoublingTime = self.parameters["nutrientToDoublingTime"]
        self.replichore_lengths = self.parameters["replichore_lengths"]
        self.sequences = self.parameters["sequences"]
        self.polymerized_dntp_weights = self.parameters["polymerized_dntp_weights"]
        self.replication_coordinate = self.parameters["replication_coordinate"]
        self.D_period = self.parameters["D_period"]
        self.replisome_protein_mass = self.parameters["replisome_protein_mass"]
        self.no_child_place_holder = self.parameters["no_child_place_holder"]
        self.basal_elongation_rate = self.parameters["basal_elongation_rate"]
        self.make_elongation_rates = self.parameters["make_elongation_rates"]

        # Sim options
        self.mechanistic_replisome = self.parameters["mechanistic_replisome"]

        # random state
        self.seed = self.parameters["seed"]
        self.random_state = np.random.RandomState(seed=self.seed)

        self.emit_unique = self.parameters.get("emit_unique", True)

        # Bulk molecule names
        self.replisome_trimers_subunits = self.parameters["replisome_trimers_subunits"]
        self.replisome_monomers_subunits = self.parameters[
            "replisome_monomers_subunits"
        ]
        self.dntps = self.parameters["dntps"]
        self.ppi = self.parameters["ppi"]

        self.ppi_idx = None

        self.debug = False

    def ports_schema(self):
        return {
            # bulk molecules
            "bulk": numpy_schema("bulk"),
            "listeners": {
                "mass": listener_schema({"cell_mass": 0.0}),
                "replication_data": listener_schema(
                    {"critical_initiation_mass": 0.0, "critical_mass_per_oriC": 0.0}
                ),
            },
            "environment": {
                "media_id": {"_default": "", "_updater": "set"},
            },
            "plasmid_active_replisomes": numpy_schema(
                "plasmid_active_replisomes", emit=self.parameters["emit_unique"]
            ),
            "oriVs": numpy_schema("oriVs", emit=self.parameters["emit_unique"]),
            "plasmid_domains": numpy_schema(
                "plasmid_domains", emit=self.parameters["emit_unique"]
            ),
            "full_plasmids": numpy_schema(
                "full_plasmids", emit=self.parameters["emit_unique"]
            ),
            "timestep": {"_default": self.parameters["time_step"]},
        }

    def calculate_request(self, timestep, states):
        if self.ppi_idx is None:
            self.ppi_idx = bulk_name_to_idx(self.ppi, states["bulk"]["id"])
            self.replisome_trimers_idx = bulk_name_to_idx(
                self.replisome_trimers_subunits, states["bulk"]["id"]
            )
            self.replisome_monomers_idx = bulk_name_to_idx(
                self.replisome_monomers_subunits, states["bulk"]["id"]
            )
            self.dntps_idx = bulk_name_to_idx(self.dntps, states["bulk"]["id"])
        requests = {}
        # Get total count of existing oriV's
        n_oriV = states["oriVs"]["_entryState"].sum()
        # If there are no origins, return immediately
        if n_oriV == 0:
            return requests

        # If replication should be initiated, request subunits required for
        # building one replisome per one origin of replication, and edit
        # access to oriC and plasmid domain attributes
        requests["bulk"] = []

        n_active_replisomes = states["plasmid_active_replisomes"]["_entryState"].sum()
        n_full_plasmids = states["full_plasmids"]["_entryState"].sum()
        # Get current locations of all replication forks
        (fork_coordinates, domain_index_replisome) = attrs(
            states["plasmid_active_replisomes"], ["coordinates", "domain_index"]
        )
        # Boolean array: True if fork is at 0 (ready to replicate)
        ready_to_replicate_mask = fork_coordinates == 0

        # At time t=1s, empty replisome attributes because we want to assemble replisome subunits
        # TODO: See whether we can make this more generalized in initial_conditions.py
        if states["global_time"] == 1:
            domain_index_replisome = []
            n_active_replisomes = 0

        # Get attributes of existing plasmid domains
        domain_index_existing_plasmid = attrs(states["full_plasmids"], ["domain_index"])
        # for newly replicated plasmids without active replisomes yet
        idle_plasmid_domains = np.setdiff1d(
            domain_index_existing_plasmid, domain_index_replisome
        )

        if len(idle_plasmid_domains) > 0:
            requests["bulk"].append(
                (self.replisome_trimers_idx, 3 * len(idle_plasmid_domains))
            )
            requests["bulk"].append(
                (self.replisome_monomers_idx, 1 * len(idle_plasmid_domains))
            )

        # If there are no active forks return

        if n_active_replisomes == 0:
            if self.debug:
                import os
                import json

                class NpEncoder(json.JSONEncoder):
                    def default(self, obj):
                        if isinstance(obj, np.integer):
                            return int(obj)
                        if isinstance(obj, np.floating):
                            return float(obj)
                        if isinstance(obj, np.ndarray):
                            return obj.tolist()
                        return super(NpEncoder, self).default(obj)

                debug_outdir = "allocator_debug_out"
                os.makedirs(debug_outdir, exist_ok=True)

                debug_record = {
                    "source": "plasmid_replication",
                    "time": states["global_time"],
                    "active replisomes": n_active_replisomes,
                    "active replisome id": domain_index_replisome,
                    "ready plasmids": domain_index_replisome[ready_to_replicate_mask],
                    "idle plasmids": idle_plasmid_domains,
                    "total plasmids": n_full_plasmids,
                }
                with open(
                    os.path.join(debug_outdir, "plasmid_allocator3.jsonl"), "a"
                ) as f:
                    f.write(json.dumps(debug_record, cls=NpEncoder) + "\n")
            return requests

        sequence_length = np.abs(np.repeat(fork_coordinates, 2))

        self.elongation_rates = self.make_elongation_rates(
            self.random_state,
            len(self.sequences),
            self.basal_elongation_rate,
            states["timestep"],
        )
        # changes made here for plasmid
        sequences = buildSequences(
            self.sequences,
            np.tile(np.arange(2), n_active_replisomes),
            sequence_length,
            self.elongation_rates,
        )

        # Count number of each dNTP in sequences for the next timestep
        sequenceComposition = np.bincount(
            sequences[sequences != polymerize.PAD_VALUE], minlength=4
        )

        # If one dNTP is limiting then limit the request for the other three by
        # the same ratio
        dNtpsTotal = counts(states["bulk"], self.dntps_idx)
        maxFractionalReactionLimit = (
            np.fmin(1, dNtpsTotal / sequenceComposition)
        ).min()

        # Request dNTPs
        requests["bulk"].append(
            (
                self.dntps_idx,
                (maxFractionalReactionLimit * sequenceComposition).astype(int),
            )
        )
        print(f"Global time: {states['global_time']}")
        if states["global_time"] >= 1327:
            print(f"Global time for debugging: {states['global_time']}")

        if self.debug:
            import os
            import json

            class NpEncoder(json.JSONEncoder):
                def default(self, obj):
                    if isinstance(obj, np.integer):
                        return int(obj)
                    if isinstance(obj, np.floating):
                        return float(obj)
                    if isinstance(obj, np.ndarray):
                        return obj.tolist()
                    return super(NpEncoder, self).default(obj)

            debug_outdir = "allocator_debug_out"
            os.makedirs(debug_outdir, exist_ok=True)

            debug_record = {
                "source": "plasmid_replication",
                "time": states["global_time"],
                "active replisomes": n_active_replisomes,
                "active replisome id": domain_index_replisome,
                "ready plasmids": domain_index_replisome[ready_to_replicate_mask],
                "idle plasmids": idle_plasmid_domains,
                "total plasmids": n_full_plasmids,
            }
            with open(os.path.join(debug_outdir, "plasmid_allocator3.jsonl"), "a") as f:
                f.write(json.dumps(debug_record, cls=NpEncoder) + "\n")

        return requests

    def evolve_state(self, timestep, states):
        # Initialize the update dictionary

        update = {
            "bulk": [],
            "plasmid_active_replisomes": {},
            "oriVs": {},
            "plasmid_domains": {},
            "full_plasmids": {},
            "listeners": {"replication_data": {}},
        }

        # Module 1: Replication initiation
        # Get number of existing replisomes and oriCs
        n_active_replisomes = states["plasmid_active_replisomes"]["_entryState"].sum()
        n_oriV = states["oriVs"]["_entryState"].sum()
        n_full_plasmids = states["full_plasmids"]["_entryState"].sum()
        # Get current fork coordinates and associated domain indices
        (fork_coordinates, domain_index_replisome) = attrs(
            states["plasmid_active_replisomes"], ["coordinates", "domain_index"]
        )

        # Boolean array: True if fork is at 0 (ready to replicate)
        # ready_to_replicate_mask = fork_coordinates == 0

        if states["global_time"] == 1:
            domain_index_replisome = []
            n_active_replisomes = 0
            # ready_to_replicate_mask = np.zeros(n_active_replisomes, dtype=bool)

        # Domain indices of plasmids ready to replicate
        # ready_domains = []
        # if np.any(ready_to_replicate_mask):
        #     ready_domains = domain_index_replisome[ready_to_replicate_mask]

        # Get attributes of existing plasmid domains
        domain_index_existing_domain, child_domains = attrs(
            states["plasmid_domains"], ["domain_index", "child_domains"]
        )
        domain_index_existing_plasmid = attrs(states["full_plasmids"], ["domain_index"])
        # for newly replicated plasmids without active replisomes yet
        idle_plasmid_domains = np.setdiff1d(
            domain_index_existing_plasmid, domain_index_replisome
        )

        # If there are no plasmids, return immediately
        if n_full_plasmids == 0:
            return update

        initiate_replication = False
        max_new_replisomes = 0
        if len(idle_plasmid_domains) > 0:
            # Get number of available replisome subunits
            n_replisome_trimers = counts(states["bulk"], self.replisome_trimers_idx)
            n_replisome_monomers = counts(states["bulk"], self.replisome_monomers_idx)

            # Calculate the maximum no.of replisomes that can be assembled in this time step
            min_trimers = int(np.min(n_replisome_trimers))
            min_monomers = int(np.min(n_replisome_monomers))
            max_by_trimers = min_trimers // 3
            max_by_monomers = min_monomers // 1

            max_new_replisomes = min(max_by_trimers, max_by_monomers)
            initiate_replication = (
                not self.mechanistic_replisome or max_new_replisomes != 0
            )

        # If all conditions are met, initiate a round of replication on max no.of orivs as possible
        if initiate_replication:
            # Get attributes of existing oriCs and domains
            (domain_index_existing_oriv,) = attrs(states["oriVs"], ["domain_index"])

            # Get indexes of the domains that would be getting child domains
            # (domains that contain an origin)
            new_parent_domains = np.where(
                np.in1d(domain_index_existing_domain, domain_index_existing_oriv)
            )[0]

            # Calculate counts of new replisomes and domains to add. changes made here for plasmid
            n_new_replisome = 0
            n_new_domain = 0
            domain_index_new = []

            if (
                not np.array_equal(idle_plasmid_domains, [0])
                and max_new_replisomes != 0
            ):
                n_new_replisome = min(len(idle_plasmid_domains), max_new_replisomes)

                n_new_domain = 2 * max_new_replisomes

                # Calculate the domain indexes of new domains and oriC's
                max_domain_index = domain_index_existing_domain.max()
                domain_index_new = np.arange(
                    max_domain_index + 1,
                    max_domain_index + 2 * max_new_replisomes + 1,
                    dtype=np.int32,
                )

            # Add new oriC's, and reset attributes of existing oriC's
            # All oriC's must be assigned new domain indexes
            if len(domain_index_new) > 0:
                if n_oriV > max_new_replisomes:
                    n_excess_orivs = n_oriV - max_new_replisomes
                    domain_index_new_oriv = np.concatenate(
                        (domain_index_existing_oriv[-n_excess_orivs:], domain_index_new)
                    )
                    update["oriVs"]["set"] = {
                        "domain_index": domain_index_new_oriv[:n_oriV]
                    }
                    update["oriVs"]["add"] = {
                        "domain_index": domain_index_new_oriv[n_oriV:],
                    }
                else:
                    update["oriVs"]["set"] = {"domain_index": domain_index_new[:n_oriV]}
                    update["oriVs"]["add"] = {
                        "domain_index": domain_index_new[n_oriV:],
                    }

            # Add and set attributes of newly created plasmid replisomes.
            # New replisomes inherit the domain indexes of the oriV's they
            # were initiated from. Only one replisome is formed per oriV (unidirectional replication).

            coordinates_replisome = np.zeros(n_new_replisome, dtype=np.int64)

            # For plasmids, replication is unidirectional, so no left/right replichore distinction.
            right_replichore = np.full(n_new_replisome, False, dtype=bool)

            # candidate_domain_index_new_replisome = np.setdiff1d(
            #     idle_plasmid_domains, domain_index_replisome
            # )

            domain_index_new_replisome = idle_plasmid_domains[:n_new_replisome]

            massDiff_protein_new_replisome = np.full(
                n_new_replisome,
                self.replisome_protein_mass if self.mechanistic_replisome else 0.0,
            )

            update["plasmid_active_replisomes"]["add"] = {
                "coordinates": coordinates_replisome,
                "right_replichore": right_replichore,
                "domain_index": domain_index_new_replisome,
                "massDiff_protein": massDiff_protein_new_replisome,
            }

            if n_new_domain != 0:
                # Add and set attributes of new plasmid domains.
                # Each new domain should have no children initially.
                new_child_domains = np.full(
                    (n_new_domain, 2), self.no_child_place_holder, dtype=np.int32
                )

                new_domains_update = {
                    "add": {
                        "domain_index": domain_index_new,
                        "child_domains": new_child_domains,
                    }
                }

                # Add new domains as children of existing domains
                # Each parent only gets one child domain (unidirectional replication)
                # child_domains[new_parent_domains, 0] = domain_index_new
                # Leave the second column (child_domains[:, 1]) as placeholder
                if new_parent_domains.size > 0:
                    if new_parent_domains.size != max_new_replisomes:
                        new_parent_domains = new_parent_domains[:max_new_replisomes]
                    child_domains[new_parent_domains] = domain_index_new.reshape(-1, 2)

                existing_domains_update = {"set": {"child_domains": child_domains}}

                update["plasmid_domains"].update(
                    {**new_domains_update, **existing_domains_update}
                )

            # Decrement counts of replisome subunits
            if self.mechanistic_replisome:
                update["bulk"].append(
                    (self.replisome_trimers_idx, -3 * n_new_replisome)
                )
                update["bulk"].append(
                    (self.replisome_monomers_idx, -1 * n_new_replisome)
                )

        # Module 2: replication elongation
        # If no active replisomes are present, return immediately
        # Note: the new replication forks added in the previous module are not
        # elongated until the next timestep.
        if n_active_replisomes == 0:
            return update

        # Get allocated counts of dNTPs
        dNtpCounts = counts(states["bulk"], self.dntps_idx)

        # Get attributes of existing replisomes
        (
            domain_index_replisome,
            right_replichore,
            coordinates_replisome,
        ) = attrs(
            states["plasmid_active_replisomes"],
            ["domain_index", "right_replichore", "coordinates"],
        )

        # Build sequences to polymerize
        sequence_length = np.abs(np.repeat(coordinates_replisome, 2))
        sequence_indexes = np.tile(np.arange(2), n_active_replisomes)

        sequences = buildSequences(
            self.sequences, sequence_indexes, sequence_length, self.elongation_rates
        )

        # Use polymerize algorithm to quickly calculate the number of
        # elongations each fork catalyzes
        reactionLimit = dNtpCounts.sum()

        active_elongation_rates = self.elongation_rates[sequence_indexes]

        result = polymerize(
            sequences,
            dNtpCounts,
            reactionLimit,
            self.random_state,
            active_elongation_rates,
        )

        sequenceElongations = result.sequenceElongation
        dNtpsUsed = result.monomerUsages

        # Compute mass increase for each elongated sequence
        mass_increase_dna = computeMassIncrease(
            sequences,
            sequenceElongations,
            self.polymerized_dntp_weights.asNumber(units.fg),
        )

        # Compute masses that should be added to each replisome
        added_dna_mass = mass_increase_dna[0::2] + mass_increase_dna[1::2]

        # Update positions of each fork
        updated_length = sequence_length + sequenceElongations
        updated_coordinates = updated_length[0::2]

        # Update attributes and submasses of replisomes
        (current_dna_mass,) = attrs(
            states["plasmid_active_replisomes"], ["massDiff_DNA"]
        )

        update["plasmid_active_replisomes"].update(
            {
                "set": {
                    "coordinates": updated_coordinates,
                    "massDiff_DNA": current_dna_mass + added_dna_mass,
                }
            }
        )

        # Update counts of polymerized metabolites
        update["bulk"].append((self.dntps_idx, -dNtpsUsed))
        update["bulk"].append((self.ppi_idx, dNtpsUsed.sum()))

        # Module 3: replication termination
        # Determine if any forks have reached the end of their sequences. If
        # so, delete the replisomes and domains that were terminated.
        # For plasmids, termination occurs when replisome reaches the end of the sequence
        terminated_replisomes = updated_coordinates >= self.replichore_lengths

        # For debugging: force termination
        # terminated_replisomes = np.ones_like(updated_coordinates, dtype=bool)

        # If any forks were terminated,
        if terminated_replisomes.sum() > 0:
            # Get domain indexes of terminated forks
            terminated_domains = np.unique(
                domain_index_replisome[terminated_replisomes]
            )

            # Get attributes of existing domains and full plasmids
            (
                domain_index_domains,
                child_domains,
            ) = attrs(states["plasmid_domains"], ["domain_index", "child_domains"])
            (domain_index_full_plasmid,) = attrs(
                states["full_plasmids"], ["domain_index"]
            )

            # Initialize array of replisomes that should be deleted
            replisomes_to_delete = np.zeros_like(domain_index_replisome, dtype=np.bool_)

            # Count number of new full plasmids that should be created
            n_new_plasmids = 0

            # Initialize array for domain indexes of new full plasmids
            domain_index_new_full_plasmid = []

            for terminated_domain_index in terminated_domains:
                # Get all terminated replisomes in the terminated domain
                terminated_domain_matching_replisomes = np.logical_and(
                    domain_index_replisome == terminated_domain_index,
                    terminated_replisomes,
                )

                # changes made here because only a single replication fork is present
                if terminated_domain_matching_replisomes.any():
                    # Tag replisomes and domains with the given domain index
                    # for deletion
                    replisomes_to_delete = np.logical_or(
                        replisomes_to_delete, terminated_domain_matching_replisomes
                    )

                    domain_mask = domain_index_domains == terminated_domain_index

                    # Get child domains of deleted domain
                    child_domains_this_domain = child_domains[
                        np.where(domain_mask)[0][0], :
                    ]

                    # Modify domain index of one existing full plasmid to
                    # index of first child domain
                    domain_index_full_plasmid = domain_index_full_plasmid.copy()
                    domain_index_full_plasmid[
                        np.where(domain_index_full_plasmid == terminated_domain_index)[
                            0
                        ]
                    ] = child_domains_this_domain[0]

                    # Increment count of new full plasmid
                    n_new_plasmids += 1

                    # Append plasmid index of new full plasmid
                    domain_index_new_full_plasmid.append(child_domains_this_domain[1])

            # Delete terminated replisomes
            update["plasmid_active_replisomes"]["delete"] = np.where(
                replisomes_to_delete
            )[0]

            # Generate new full plasmid molecules
            if n_new_plasmids > 0:
                plasmid_add_update = {
                    "add": {
                        "domain_index": domain_index_new_full_plasmid,
                        "division_time": states["global_time"] + self.D_period,
                        "has_triggered_division": False,
                    }
                }

                # Reset domain index of existing plasmids that have finished
                # replication
                plasmid_existing_update = {
                    "set": {"domain_index": domain_index_full_plasmid}
                }

                update["full_plasmids"].update(
                    {**plasmid_add_update, **plasmid_existing_update}
                )

            # Increment counts of replisome subunits
            if self.mechanistic_replisome:
                update["bulk"].append(
                    (self.replisome_trimers_idx, 3 * replisomes_to_delete.sum())
                )
                update["bulk"].append(
                    (self.replisome_monomers_idx, replisomes_to_delete.sum())
                )

        return update


def test_plasmid_replication():
    from ecoli.library.sim_data import LoadSimData

    sim_data_default = "../../out/plasmidwithsequence/parca/kb/simData.cPickle"
    load_sim_data = LoadSimData(sim_data_default)

    replication_config = load_sim_data.get_plasmid_replication_config()

    # the full initial state
    initial_state = load_sim_data.generate_initial_state()

    # test_config = {}
    process = PlasmidReplication(replication_config)
    # assert process is not None

    # interval = 1
    input_state1 = {
        "bulk": initial_state["bulk"],
        "environment": initial_state["environment"],
        "plasmid_active_replisomes": initial_state["unique"][
            "plasmid_active_replisome"
        ],
        "oriVs": initial_state["unique"]["oriV"],
        "plasmid_domains": initial_state["unique"]["plasmid_domain"],
        "full_plasmids": initial_state["unique"]["full_plasmid"],
        "listeners": {"mass": {"cell_mass": 2000.0}},  # arbitrary fg to start,
        "timestep": replication_config["time_step"],
    }

    requests_history = []
    bulk_ids = input_state1["bulk"]["id"]

    import pandas as pd

    # Call calculate_request
    requests = process.calculate_request(input_state1["timestep"], input_state1)

    # Convert requests["bulk"] into a flat dictionary
    bulk_dict = {}
    for ids, count in requests.get("bulk", []):
        ids = np.atleast_1d(ids)
        if np.isscalar(count):
            count = np.full_like(ids, count, dtype=int)
        else:
            count = np.array(count, dtype=int)
        for idx, cnt in zip(ids, count):
            mol_name = bulk_ids[idx]  # map from id to molecule name
            bulk_dict[mol_name] = int(cnt)

    # Append to history
    requests_history.append({**bulk_dict})

    # Convert to dataframe
    df_requests = pd.DataFrame(requests_history)

    # plotting
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Melt the dataframe to long format for easier plotting
    df_long = df_requests.melt(var_name="Molecule", value_name="Count")

    plt.figure(figsize=(10, 6))
    sns.barplot(data=df_long, x="Molecule", y="Count", palette="viridis")

    plt.xlabel("Bulk Molecule")
    plt.ylabel("Requested Count")
    plt.title("Bulk Molecule Requests for Plasmid at current time step")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig("plasmid_bulk_requests_plot.png", dpi=300)

    # call evolve_state
    update = process.evolve_state(input_state1["timestep"], input_state1)
    print(update.keys())


if __name__ == "__main__":
    test_plasmid_replication()
