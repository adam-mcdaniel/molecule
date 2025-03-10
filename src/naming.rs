use super::*;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use std::collections::{BTreeMap, HashMap, HashSet};
use tracing::*; // added logging from tracing

// --- NEW: Global functional group detection functions ---

/// Detect an ester functional group anywhere in the molecule.
/// Returns (acyl_carbon, alkoxy_oxygen) if found.
fn detect_ester_global(graph: &MoleculeGraph) -> Option<(NodeIndex, NodeIndex)> {
    debug!("Running detect_ester_global");
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            let mut has_carbonyl = false;
            let mut alkoxy_oxygen: Option<NodeIndex> = None;
            for nbr in graph.neighbors(node) {
                if graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        match graph.edge_weight(edge) {
                            Some(Bond::Double) => {
                                has_carbonyl = true;
                                debug!("Found carbonyl on node {:?} via neighbor {:?}", node, nbr);
                            }
                            Some(Bond::Single) => {
                                // Check if oxygen connects to a carbon (other than node)
                                for nbr2 in graph.neighbors(nbr) {
                                    if nbr2 != node && (graph[nbr2].is_carbon()) {
                                        alkoxy_oxygen = Some(nbr);
                                        debug!("Found alkoxy oxygen {:?} connected to node {:?} from {:?}", nbr, node, nbr2);
                                        break;
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
            if has_carbonyl && alkoxy_oxygen.is_some() {
                debug!(
                    "Ester detected: acyl carbon {:?} with alkoxy oxygen {:?}",
                    node,
                    alkoxy_oxygen.unwrap()
                );
                return Some((node, alkoxy_oxygen.unwrap()));
            }
        }
    }
    debug!("No ester group found in detect_ester_global");
    None
}

/// Detect an aldehyde group (terminal carbonyl) in the molecule.
/// Returns the longest carbon chain and the index (in that chain) of the aldehyde carbon.
fn detect_aldehyde_global(graph: &MoleculeGraph) -> Option<(Vec<NodeIndex>, usize)> {
    debug!("Running detect_aldehyde_global");
    let main_chain = find_longest_carbon_chain(graph);
    if main_chain.is_empty() {
        debug!("No carbon chain found for aldehyde detection");
        return None;
    }
    // Check first carbon
    let first = main_chain[0];
    for nbr in graph.neighbors(first) {
        if !main_chain.contains(&nbr) && graph[nbr] == Element::O {
            if let Some(edge) = graph.find_edge(first, nbr) {
                if graph.edge_weight(edge) == Some(&Bond::Double) {
                    debug!("Aldehyde detected at start of chain: node {:?}", first);
                    return Some((main_chain, 0));
                }
            }
        }
    }
    // Check last carbon
    let last_index = main_chain.len() - 1;
    let last = main_chain[last_index];
    for nbr in graph.neighbors(last) {
        if !main_chain.contains(&nbr) && graph[nbr] == Element::O {
            if let Some(edge) = graph.find_edge(last, nbr) {
                if graph.edge_weight(edge) == Some(&Bond::Double) {
                    debug!("Aldehyde detected at end of chain: node {:?}", last);
                    return Some((main_chain, last_index));
                }
            }
        }
    }
    debug!("No aldehyde group found in detect_aldehyde_global");
    None
}

/// Helper: performs DFS (only over carbon atoms) starting from `start`,
/// ignoring the specified `exclude` node (typically the oxygen).
fn dfs_component(graph: &MoleculeGraph, start: NodeIndex, exclude: NodeIndex) -> Vec<NodeIndex> {
    let mut visited = HashSet::new();
    let mut stack = vec![start];
    while let Some(n) = stack.pop() {
        if visited.insert(n) {
            for nbr in graph.neighbors(n) {
                if nbr == exclude {
                    continue;
                }
                if graph[nbr].is_carbon() {
                    if !visited.contains(&nbr) {
                        stack.push(nbr);
                    }
                }
            }
        }
    }
    visited.into_iter().collect()
}

/// Detect an ether group in the molecule (i.e. a single oxygen not in a carbonyl
/// that bridges two carbon fragments). If found, returns a tuple:
/// (main_chain, alkoxy_chain, oxygen) where the main_chain is chosen as the larger fragment.
fn detect_ether_global(
    graph: &MoleculeGraph,
) -> Option<(Vec<NodeIndex>, Vec<NodeIndex>, NodeIndex)> {
    debug!("Running detect_ether_global");
    // Look for oxygen nodes in the graph.
    for o in graph.node_indices().filter(|&n| graph[n] == Element::O) {
        let mut is_carbonyl = false;
        let mut carbon_neighbors = Vec::new();
        for nbr in graph.neighbors(o) {
            if graph[nbr].is_carbon() {
                if let Some(edge) = graph.find_edge(o, nbr) {
                    if graph.edge_weight(edge) == Some(&Bond::Double) {
                        is_carbonyl = true;
                        debug!(
                            "Oxygen node {:?} is part of a carbonyl with neighbor {:?}",
                            o, nbr
                        );
                        break;
                    } else {
                        carbon_neighbors.push(nbr);
                    }
                }
            }
        }
        if is_carbonyl {
            continue;
        }
        if carbon_neighbors.len() == 2 {
            let comp1 = dfs_component(graph, carbon_neighbors[0], o);
            let comp2 = dfs_component(graph, carbon_neighbors[1], o);
            debug!(
                "Ether candidate at oxygen {:?}: comp1 size {}, comp2 size {}",
                o,
                comp1.len(),
                comp2.len()
            );
            if comp1.len() >= comp2.len() {
                return Some((comp1, comp2, o));
            } else {
                return Some((comp2, comp1, o));
            }
        }
    }
    debug!("No ether group found in detect_ether_global");
    None
}

// --- NEW: Exocyclic Aldehyde Detection ---
/// If a carbonyl (aldehyde) carbon is exocyclic—that is, not part of a ring but attached
/// to a heterocycle—this function returns its proper name using the heterocyclic acid rules.
fn detect_exocyclic_aldehyde(graph: &MoleculeGraph) -> Option<String> {
    debug!("Running detect_exocyclic_aldehyde");
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            for nbr in graph.neighbors(node) {
                if graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        if graph.edge_weight(edge) == Some(&Bond::Double) {
                            // node is a potential carbonyl carbon.
                            for nbr2 in graph.neighbors(node) {
                                if nbr2 != nbr {
                                    if let Some(ring) = find_ring_containing_node(graph, nbr2) {
                                        trace!("Found ring containing node while searching for exocyclic aldehyde {:?}", nbr2);
                                        if let Some(name) = name_exocyclic_aldehyde(graph, node) {
                                            debug!("Exocyclic aldehyde found at node {:?} attached to ring", node);
                                            return Some(name);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    debug!("No exocyclic aldehyde detected");
    None
}

/// Names a single ester molecule using the globally detected ester atoms.
fn iupac_ester_name(
    graph: &MoleculeGraph,
    acyl_carbon: NodeIndex,
    alkoxy_oxygen: NodeIndex,
) -> String {
    debug!(
        "Naming ester: acyl carbon {:?}, alkoxy oxygen {:?}",
        acyl_carbon, alkoxy_oxygen
    );
    let acyl_name = name_acyl_group(graph, acyl_carbon);
    if let Some(alkoxy_carb) = graph
        .neighbors(alkoxy_oxygen)
        .find(|&nbr| nbr != acyl_carbon && graph[nbr].is_carbon())
    {
        if let Some(ring_name) = alkoxy_ring_name(graph, alkoxy_carb) {
            debug!("Using aromatic alkoxy ring name: {}", ring_name);
            return format!("{}-{}", ring_name, acyl_name);
        }
        if is_tert_butyl(graph, alkoxy_carb, alkoxy_oxygen) {
            debug!("Detected tert-butyl substituent");
            return format!("tert-butyl-{}", acyl_name);
        }
        // NEW: Try iso substituent naming first.
        if let Some(iso) = detect_iso_substituent(graph, alkoxy_carb, alkoxy_oxygen) {
            debug!("Detected iso substituent: {}", iso);
            return format!("{}-{}", iso, acyl_name);
        }
        // Fallback: use DFS branch size.
        let branch_size = count_branch_carbons(graph, alkoxy_carb, alkoxy_oxygen);
        let alkoxy_name = match branch_size {
            1 => "methyl",
            2 => "ethyl",
            3 => "propyl",
            4 => "butyl",
            5 => "pentyl",
            _ => "alkyl",
        };
        debug!("Ester naming: using fallback alkoxy prefix {}", alkoxy_name);
        return format!("{}-{}", alkoxy_name, acyl_name);
    }
    acyl_name
}

/// Names an aldehyde molecule given its longest carbon chain and the index
/// (in that chain) at which the terminal aldehyde group occurs.
fn iupac_aldehyde_name(
    graph: &MoleculeGraph,
    mut main_chain: Vec<NodeIndex>,
    aldehyde_index: usize,
) -> String {
    debug!(
        "Naming aldehyde: chain length {}, aldehyde index {}",
        main_chain.len(),
        aldehyde_index
    );
    if let Some(exo_aldehyde) = name_exocyclic_aldehyde(graph, main_chain[aldehyde_index]) {
        debug!("Using exocyclic aldehyde naming: {}", exo_aldehyde);
        return exo_aldehyde;
    }
    if aldehyde_index != 0 {
        main_chain.reverse();
        debug!("Reversed main chain to place aldehyde at start");
    }
    let mut subs = identify_substituents(graph, &main_chain);
    subs.retain(|(pos, name)| !(*pos == 1 && name == "oxo"));
    let numbering_order = choose_numbering(&main_chain, &subs);
    let final_chain = match numbering_order {
        NumberingOrder::Original => main_chain.clone(),
        NumberingOrder::Reversed => {
            let mut rev = main_chain.clone();
            rev.reverse();
            rev
        }
    };
    let subs_final = assign_substituent_positions(&final_chain, &subs, numbering_order);
    let substituent_prefix = format_substituents(&subs_final);
    let chain_length = final_chain.len();
    let root = match chain_length {
        1 => "methan",
        2 => "ethan",
        3 => "propan",
        4 => "butan",
        5 => "pentan",
        6 => "hexan",
        7 => "heptan",
        8 => "octan",
        9 => "nonan",
        10 => "decan",
        n => return format!("{}-carbon chain (not fully supported)", n),
    };
    debug!(
        "Aldehyde naming: final chain length {}, substituent prefix '{}'",
        chain_length, substituent_prefix
    );
    format!("{}al", substituent_prefix + root)
}

/// Names an ether molecule given the two disconnected carbon fragments and the bridging oxygen.
/// The main_chain is chosen as the larger (or preferred) fragment.
fn iupac_ether_name(
    graph: &MoleculeGraph,
    mut main_chain: Vec<NodeIndex>,
    alkoxy_chain: Vec<NodeIndex>,
    oxygen: NodeIndex,
) -> String {
    debug!(
        "Naming ether: main chain length {}, alkoxy chain length {}, oxygen {:?}",
        main_chain.len(),
        alkoxy_chain.len(),
        oxygen
    );
    if let Some(&attach) = main_chain
        .iter()
        .find(|&&n| graph.find_edge(n, oxygen).is_some())
    {
        while main_chain[0] != attach {
            main_chain.rotate_left(1);
        }
        let main_name = match main_chain.len() {
            1 => "methane",
            2 => "ethane",
            3 => "propane",
            4 => "butane",
            5 => "pentane",
            6 => "hexane",
            7 => "heptane",
            8 => "octane",
            9 => "nonane",
            10 => "decane",
            n => return format!("{}-carbon chain (not fully supported)", n),
        };
        let ether_prefix = match alkoxy_chain.len() {
            1 => "methoxy",
            2 => "ethoxy",
            3 => "propoxy",
            4 => "butoxy",
            5 => "pentoxy",
            _ => "alkoxy",
        };
        debug!("Ether naming: using prefix '{}'", ether_prefix);
        return format!("1-{}-{}", ether_prefix, main_name);
    }
    debug!("Fallback to standard acyclic naming for ether");
    iupac_acyclic_name(graph)
}

/// Detect all ester functional groups in the molecule.
/// Returns a vector of (acyl_carbon, alkoxy_oxygen) pairs.
fn detect_esters_global(graph: &MoleculeGraph) -> Vec<(NodeIndex, NodeIndex)> {
    let mut esters = Vec::new();
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            let mut has_carbonyl = false;
            let mut alkoxy_oxygens = Vec::new();
            for nbr in graph.neighbors(node) {
                if graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        match graph.edge_weight(edge) {
                            Some(Bond::Double) => {
                                has_carbonyl = true;
                            }
                            Some(Bond::Single) => {
                                // Check if oxygen connects to a carbon (other than node)
                                for nbr2 in graph.neighbors(nbr) {
                                    if nbr2 != node && graph[nbr2].is_carbon() {
                                        alkoxy_oxygens.push(nbr);
                                        break;
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
            if has_carbonyl && !alkoxy_oxygens.is_empty() {
                for alkoxy in alkoxy_oxygens {
                    esters.push((node, alkoxy));
                }
            }
        }
    }
    esters
}

/// Modified general acyl group naming functi1on for esters.
fn name_acyl_group(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> String {
    // First, if the acyl carbon is directly bonded to an aromatic atom that is part
    // of a six-membered aromatic ring, assume benzoic acid.
    let is_benzo = graph.neighbors(acyl_carbon).any(|nbr| {
        if graph[nbr].is_carbon() {
            if let Some(ring) = find_ring(graph) {
                ring.contains(&nbr) && ring.len() == 6 && ring.iter().all(|&n| graph[n].is_carbon())
            } else {
                false
            }
        } else {
            false
        }
    });
    if is_benzo {
        return "benzoate".to_string();
    }

    // Next, check if the acyl carbon is part of a ring.
    if let Some(ring) = find_ring_containing_node(graph, acyl_carbon) {
        trace!("Found ring containing acyl carbon: {:?}", ring);
        return name_heterocyclic_acid(graph, &ring, acyl_carbon);
    }

    // Otherwise, check if one of its non‑oxygen neighbors is in a heterocyclic ring.
    if let Some(exo_name) = name_exocyclic_acid(graph, acyl_carbon) {
        return exo_name;
    }
    if let Some(name) = name_exocyclic_aldehyde(graph, acyl_carbon) {
        return name;
    }

    // Fallback: treat the acyl group as acyclic.
    let acyl_chain = find_acyl_chain(graph, acyl_carbon);
    // --- NEW: Special-case for oxalic acid ---
    if acyl_chain.len() == 2 {
        let mut count = 0;
        for &node in &acyl_chain {
            for nbr in graph.neighbors(node) {
                if !acyl_chain.contains(&nbr) && graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        if graph.edge_weight(edge) == Some(&Bond::Double) {
                            count += 1;
                            break;
                        }
                    }
                }
            }
        }
        if count == 2 {
            return "oxalate".to_string();
        }
    }
    // NEW: Incorporate substituent naming into the acid fragment.
    // Identify substituents on the acyl chain...
    let subs = identify_substituents(graph, &acyl_chain);
    // ...and remove any substituent at position 1 (the carbonyl carbon)
    let subs: Vec<_> = subs.into_iter().filter(|(pos, _)| *pos != 1).collect();

    let numbering_order = choose_numbering(&acyl_chain, &subs);
    let final_acyl_chain = match numbering_order {
        NumberingOrder::Original => acyl_chain.clone(),
        NumberingOrder::Reversed => {
            let mut rev = acyl_chain.clone();
            rev.reverse();
            rev
        }
    };
    let subs_final = assign_substituent_positions(&final_acyl_chain, &subs, numbering_order);
    let substituent_prefix = format_substituents(&subs_final);
    let chain_length = final_acyl_chain.len();
    let root = match chain_length {
        1 => "methan",
        2 => "ethan",
        3 => "propan",
        4 => "butan",
        5 => "pentan",
        6 => "hexan",
        7 => "heptan",
        8 => "octan",
        9 => "nonan",
        10 => "decan",
        n => return format!("{}-carbon acid chain (not fully supported)", n),
    };
    // format!("{}oate", substituent_prefix + root)
    format!("{}oate", format!("{}{}", substituent_prefix.trim_end_matches('-'), root))
}

/// Names a molecule with multiple ester groups (polyester), for example producing
/// "diethyl-propanedioate" when there are two identical ester groups.
fn iupac_polyester_name(graph: &MoleculeGraph, esters: &[(NodeIndex, NodeIndex)]) -> String {
    // For simplicity, assume all ester groups belong to the same acid fragment.
    let mut acyl_set = std::collections::HashSet::new();
    let mut alkoxy_names = Vec::new();
    for &(acyl, alkoxy) in esters {
        acyl_set.insert(acyl);
        if let Some(alkoxy_carb) = graph
            .neighbors(alkoxy)
            .find(|&nbr| nbr != acyl && graph[nbr].is_carbon())
        {
            if is_tert_butyl(graph, alkoxy_carb, alkoxy) {
                alkoxy_names.push("tert-butyl".to_string());
            } else {
                let branch_size = count_branch_carbons(graph, alkoxy_carb, alkoxy);
                let name = match branch_size {
                    1 => "methyl".to_string(),
                    2 => "ethyl".to_string(),
                    3 => "propyl".to_string(),
                    4 => "butyl".to_string(),
                    5 => "pentyl".to_string(),
                    _ => "alkyl".to_string(),
                };
                alkoxy_names.push(name);
            }
        }
    }
    // Determine how many ester groups we have.
    let ester_count = esters.len();
    // If all alkoxy groups are identical, use a multiplier prefix.
    let first_name = &alkoxy_names[0];
    let all_same = alkoxy_names.iter().all(|name| name == first_name);
    let alkoxy_prefix = if all_same {
        let multiplier = match ester_count {
            2 => "di",
            3 => "tri",
            4 => "tetra",
            5 => "penta",
            _ => "",
        };
        format!("{}{}", multiplier, first_name)
    } else {
        // If not identical, join them with commas.
        alkoxy_names.join(",")
    };
    // Use one of the acyl carbons to define the acid fragment.
    let acyl_carbon = *acyl_set.iter().next().unwrap();
    let acyl_name = name_acyl_group(graph, acyl_carbon);
    // --- NEW: Preserve "oxalate" as is ---
    let acid_name = if acyl_name == "oxalate" {
        acyl_name
    } else if acyl_name.ends_with("oate") {
        let base = &acyl_name[..acyl_name.len() - 4];
        let multiplier_suffix = match ester_count {
            2 => "dioate",
            3 => "trioate",
            4 => "tetraoate",
            5 => "pentaoate",
            _ => "oate",
        };
        format!("{}{}", base, multiplier_suffix)
    } else {
        acyl_name
    };
    format!("{}-{}", alkoxy_prefix, acid_name)
}

/// Returns the longest chain of carbons that are NOT part of any ring.
/// If none is found, returns an empty vector.
fn find_longest_acyclic_carbon_chain(graph: &MoleculeGraph) -> Vec<NodeIndex> {
    let mut longest = Vec::new();
    for node in graph.node_indices() {
        if graph[node] == Element::C && find_ring_containing_node(graph, node).is_none() {
            let mut visited = HashSet::new();
            // Use your DFS routine but then filter out any nodes that turn out to be in a ring.
            let chain = dfs_longest_path(graph, node, &mut visited)
                .into_iter()
                .filter(|&n| find_ring_containing_node(graph, n).is_none())
                .collect::<Vec<_>>();
            if chain.len() > longest.len() {
                longest = chain;
            }
        }
    }
    longest
}

/// --- Modified acyclic naming routine ---
pub fn iupac_acyclic_name(graph: &MoleculeGraph) -> String {
    debug!("Starting acyclic naming routine");
    let esters = detect_esters_global(graph);
    if !esters.is_empty() {
        debug!("Acyclic naming: Ester branch taken");
        if esters.len() == 1 {
            let (acyl_carbon, alkoxy_oxygen) = esters[0];
            return iupac_single_ester_name(graph, acyl_carbon, alkoxy_oxygen);
        } else {
            return iupac_polyester_name(graph, &esters);
        }
    }
    // NEW: Check if a free acid group is present.
    if let Some(_) = detect_acid_global(graph) {
        return iupac_acid_name(graph);
    }
    if let Some((main_chain, aldehyde_index)) = detect_aldehyde_global(graph) {
        debug!("Acyclic naming: Aldehyde branch taken");
        return iupac_aldehyde_name(graph, main_chain, aldehyde_index);
    }
    if let Some((main_chain, alkoxy_chain, oxygen)) = detect_ether_global(graph) {
        debug!("Acyclic naming: Ether branch taken");
        return iupac_ether_name(graph, main_chain, alkoxy_chain, oxygen);
    }
    if let Some(name) = detect_exocyclic_aldehyde(graph) {
        debug!("Acyclic naming: Exocyclic aldehyde branch taken");
        return name;
    }

    debug!("Acyclic naming: Falling back to standard acyclic naming");
    // Prefer an acyclic carbon chain if available.
    let main_chain = {
        let acyclic = find_longest_acyclic_carbon_chain(graph);
        if !acyclic.is_empty() { acyclic } else { find_longest_carbon_chain(graph) }
    };
    if main_chain.is_empty() {
        debug!("No valid carbon chain found; returning Unknown molecule");
        return "Unknown molecule".to_string();
    }

    let subs_original = identify_substituents(graph, &main_chain);
    let numbering_order = choose_numbering(&main_chain, &subs_original);
    let final_main_chain = match numbering_order {
        NumberingOrder::Original => main_chain.clone(),
        NumberingOrder::Reversed => {
            let mut rev = main_chain.clone();
            rev.reverse();
            rev
        }
    };
    let subs_final =
        assign_substituent_positions(&final_main_chain, &subs_original, numbering_order);
    let substituent_prefix = format_substituents(&subs_final);

    let chain_length = final_main_chain.len();
    let base_name = match chain_length {
        1 => "meth",
        2 => "eth",
        3 => "prop",
        4 => "but",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        n => return format!("{}-carbon chain (not fully supported)", n),
    };

    let mut unsat_locs = Vec::new();
    let mut double_count = 0;
    let mut triple_count = 0;
    for i in 0..(final_main_chain.len() - 1) {
        let n1 = final_main_chain[i];
        let n2 = final_main_chain[i + 1];
        if let Some(edge) = graph.find_edge(n1, n2) {
            match graph.edge_weight(edge) {
                Some(Bond::Double) => {
                    if chain_length > 2 {
                        unsat_locs.push(i + 1);
                    }
                    double_count += 1;
                }
                Some(Bond::Triple) => {
                    if chain_length > 2 {
                        unsat_locs.push(i + 1);
                    }
                    triple_count += 1;
                }
                _ => {}
            }
        }
    }
    let suffix = if double_count > 0 {
        match double_count {
            1 => "ene",
            2 => "diene",
            3 => "triene",
            _ => "ene",
        }
    } else if triple_count > 0 {
        match triple_count {
            1 => "yne",
            2 => "diyne",
            3 => "triyne",
            _ => "yne",
        }
    } else {
        "ane"
    };
    let unsaturated_base = if (double_count > 1 || triple_count > 1) && chain_length >= 4 {
        format!("{}a", base_name)
    } else {
        base_name.to_string()
    };
    let unsat_prefix = if chain_length == 2 || (double_count == 0 && triple_count == 0) {
        "".to_string()
    } else {
        let locants: Vec<String> = unsat_locs.iter().map(|loc| loc.to_string()).collect();
        format!("{}-", locants.join(","))
    };
    debug!(
        "Acyclic naming: final name composed with prefix '{}' and suffix '{}'",
        substituent_prefix, suffix
    );
    format!(
        "{}{}{}{}",
        substituent_prefix, unsat_prefix, unsaturated_base, suffix
    )
}

// --- NEW: Generic ring detection that does not filter out heteroatoms ---
fn find_ring(graph: &MoleculeGraph) -> Option<Vec<NodeIndex>> {
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    for node in graph.node_indices() {
        if !visited.contains(&node) {
            if let Some(cycle) = dfs_find_cycle_generic(graph, None, node, None, &mut stack, &mut visited)
            {
                return Some(cycle);
            }
        }
    }
    None
}

fn dfs_find_cycle_generic(
    graph: &MoleculeGraph,
    start: Option<NodeIndex>,
    current: NodeIndex,
    parent: Option<NodeIndex>,
    stack: &mut Vec<NodeIndex>,
    visited: &mut HashSet<NodeIndex>,
) -> Option<Vec<NodeIndex>> {
    visited.insert(current);
    stack.push(current);
    trace!("Visiting node: {:?}", current);
    for neighbor in graph.neighbors(current) {
        if Some(neighbor) == parent {
            continue;
        }

        if stack.contains(&neighbor) {
            match start {
                // Confirm the cycle contains the start node.
                Some(start) => {
                    if stack.contains(&start) {
                        trace!("Cycle detected with start node: {:?}", start);
                    } else {
                        continue;
                    }
                }
                None => {
                    trace!("Cycle detected without start node");
                }
            }
            // info!("Cycle detected at node: {:?}", neighbor);
            let pos = stack.iter().position(|&x| x == neighbor).unwrap();
            return Some(stack[pos..].to_vec());
        } else if !visited.contains(&neighbor) {
            if let Some(cycle) =
                dfs_find_cycle_generic(graph, start, neighbor, Some(current), stack, visited)
            {
                return Some(cycle);
            }
        }
    }
    stack.pop();
    None
}

/// Modified heterocyclic naming function.
/// If the ring is pyrrolidine, return the canonical SMILES.
fn iupac_heterocyclo_name(graph: &MoleculeGraph, ring: &Vec<NodeIndex>) -> String {
    // NEW: Check if the ring is a saturated 5-membered ring with 4 carbons and 1 nitrogen.
    if ring.len() == 5 {
        let mut count_c = 0;
        let mut count_n = 0;
        for &n in ring {
            match graph[n] {
                Element::C => count_c += 1,
                Element::N => count_n += 1,
                _ => {}
            }
        }
        if count_c == 4 && count_n == 1 {
            return "pyrrolidine".to_string();
        }
    }

    // Existing code for heterocyclic naming:
    let n = ring.len();
    let base_name = match n {
        3 => "cyclopropane",
        4 => "cyclobutane",
        5 => "cyclopentane",
        6 => "cyclohexane",
        7 => "cycloheptane",
        8 => "cyclooctane",
        9 => "cyclononane",
        10 => "cyclodecane",
        _ => return format!("{}-membered ring (not supported)", n),
    };

    let (_best_order, hetero_atoms) = best_heterocycle_numbering(ring, graph);
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (loc, elem) in hetero_atoms {
        let group_name = match elem {
            Element::O => "oxa".to_string(),
            Element::N => "aza".to_string(),
            Element::S => "thia".to_string(),
            _ => "hetero".to_string(),
        };
        groups.entry(group_name).or_default().push(loc);
    }
    let prefix = if groups.len() == 1 {
        let (group, mut locants) = groups.into_iter().next().unwrap();
        locants.sort();
        let locant_str = locants
            .iter()
            .map(|l| l.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let count = locants.len();
        let multiplier = match count {
            1 => "",
            2 => "di",
            3 => "tri",
            4 => "tetra",
            _ => "",
        };
        format!("{}-{}{}", locant_str, multiplier, group)
    } else {
        String::new()
    };

    format!("{}cyclo{}", prefix, &base_name[5..])
}

// --- NEW: Numbering for heterocyclic rings ---
fn best_heterocycle_numbering(
    ring: &[NodeIndex],
    graph: &MoleculeGraph,
) -> (Vec<NodeIndex>, Vec<(usize, Element)>) {
    let n = ring.len();
    let mut best_order = ring.to_vec();
    let mut best_locants: Vec<usize> = ring
        .iter()
        .enumerate()
        .filter(|&(_, &node)| !graph[node].is_carbon())
        .map(|(i, _)| i + 1)
        .collect();
    best_locants.sort();

    let mut best_hetero: Vec<(usize, Element)> = ring
        .iter()
        .enumerate()
        .filter(|&(_, &node)| !graph[node].is_carbon())
        .map(|(i, &node)| (i + 1, graph[node]))
        .collect();
    best_hetero.sort_by_key(|&(loc, _)| loc);

    // Try all rotations and directions.
    for start in 0..n {
        let candidate: Vec<NodeIndex> = ring[start..]
            .iter()
            .chain(ring[..start].iter())
            .cloned()
            .collect();
        for &reverse in &[false, true] {
            let candidate_order = if reverse {
                let mut rev = candidate.clone();
                rev.reverse();
                rev
            } else {
                candidate.clone()
            };
            let mut candidate_locs: Vec<usize> = Vec::new();
            let mut candidate_hetero: Vec<(usize, Element)> = Vec::new();
            for (i, &node) in candidate_order.iter().enumerate() {
                if !graph[node].is_carbon() {
                    candidate_locs.push(i + 1);
                    candidate_hetero.push((i + 1, graph[node]));
                }
            }
            candidate_locs.sort();
            if candidate_locs < best_locants {
                best_locants = candidate_locs;
                best_order = candidate_order;
                best_hetero = candidate_hetero;
            }
        }
    }
    (best_order, best_hetero)
}

#[derive(Debug, Clone, Copy)]
enum NumberingOrder {
    Original,
    Reversed,
}

struct Substituent {
    position: usize, // 1-indexed position on the parent chain.
    name: String,
}

fn find_longest_carbon_chain(graph: &MoleculeGraph) -> Vec<NodeIndex> {
    let mut longest_chain = Vec::new();
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            let mut visited = HashSet::new();
            let chain = dfs_longest_path(graph, node, &mut visited);
            if chain.len() > longest_chain.len() {
                longest_chain = chain;
            }
        }
    }
    longest_chain
}

/// For each carbon in the main chain, identify any neighbor (not in the chain)
/// as a substituent.
fn identify_substituents(graph: &MoleculeGraph, main_chain: &[NodeIndex]) -> Vec<(usize, String)> {
    let mut subs = Vec::new();
    let main_chain_set: HashSet<_> = main_chain.iter().cloned().collect();
    for (i, &chain_node) in main_chain.iter().enumerate() {
        for neighbor in graph.neighbors(chain_node) {
            if !main_chain_set.contains(&neighbor) {
                let sub_name = name_substituent(graph, neighbor, chain_node);
                subs.push((i + 1, sub_name));
            }
        }
    }
    subs
}

/// For each node in the ring, identify neighbors not in the ring as substituents.
fn identify_substituents_on_ring(
    graph: &MoleculeGraph,
    ring: &Vec<NodeIndex>,
) -> Vec<(usize, String)> {
    let mut subs = Vec::new();
    let ring_set: HashSet<_> = ring.iter().cloned().collect();
    for (i, &ring_node) in ring.iter().enumerate() {
        for neighbor in graph.neighbors(ring_node) {
            if !ring_set.contains(&neighbor) {
                let sub_name = name_substituent(graph, neighbor, ring_node);
                subs.push((i + 1, sub_name)); // positions are 1-indexed on the ring
            }
        }
    }
    subs
}

// NEW: Add a function to detect guanidino groups as substituents.
fn detect_guanidino_substituent(graph: &MoleculeGraph, start: NodeIndex, exclude: NodeIndex) -> Option<String> {
    // Look for a carbon neighbor of 'start' (excluding the parent 'exclude').
    let central_carbon = graph.neighbors(start)
        .find(|&n| n != exclude && graph[n] == Element::C)?;
    // Count the number of nitrogen neighbors (other than 'start') attached to this central carbon.
    let n_count = graph.neighbors(central_carbon)
        .filter(|&n| n != start && graph[n] == Element::N)
        .count();
    if n_count >= 2 {
        return Some("guanidino".to_string());
    }
    None
}

// NEW: Add a helper to detect a sulfhydryl substituent.
fn detect_sulfhydryl_substituent(
    graph: &MoleculeGraph,
    start: NodeIndex,
    exclude: NodeIndex,
) -> Option<String> {
    // Look for a sulfur neighbor of the substituent carbon 'start'
    // (ignoring the parent 'exclude').
    let sulfur_neighbors: Vec<_> = graph
        .neighbors(start)
        .filter(|&n| n != exclude && graph[n] == Element::S)
        .collect();
    if sulfur_neighbors.len() == 1 {
        let s = sulfur_neighbors[0];
        // Optionally: check that S has no other carbon neighbors (i.e. it is terminal)
        let other_neighbors: Vec<_> = graph
            .neighbors(s)
            .filter(|&n| n != start)
            .collect();
        if other_neighbors.iter().all(|&n| !graph[n].is_carbon()) {
            return Some("sulfhydryl".to_string());
        }
    }
    None
}

/// Returns a base name for a ring based on its elemental composition.
/// Extended to recognize pyrrole (aromatic) and pyrrolidine (saturated) for 5‑membered rings.
fn base_name_for_ring(graph: &MoleculeGraph, ring: &[NodeIndex]) -> Option<String> {
    let mut count_c = 0;
    let mut count_o = 0;
    let mut count_n = 0;
    let mut count_s = 0;
    for &n in ring {
        let elem = graph[n];
        if elem.is_carbon() {
            count_c += 1;
        } else if elem.is_oxygen() {
            count_o += 1;
        } else if elem.is_nitrogen() {
            count_n += 1;
        } else if elem.is_sulfur() {
            count_s += 1;
        }
    }
    // Recognize benzene: 6-membered ring, 6 carbons.
    if ring.len() == 6 && count_c == 6 && count_o == 0 && count_n == 0 && count_s == 0 {
        // Check if all the bonds are single. If so, it's cyclohexane.
        if ring.iter().all(|&n| graph[n].is_aromatic()) {
            return Some("benzene".to_string());
        }
        // If all the bonds are single
        if ring.iter().all(|&n| {
            // Get the edges connected to the node
            let edges = graph.edges(n);
            
            // Check if all edges are single bonds
            // Filter only edges that are to the ring
            edges.filter(|edge| ring.contains(&edge.source()) && ring.contains(&edge.target()))
                .all(|edge| edge.weight() == &Bond::Single)
        }) {
            return Some("cyclohexane".to_string());
        }

        return Some("benzene".to_string());
    }
    // Recognize pyridine: 6-membered ring, 5 carbons and 1 nitrogen.
    if ring.len() == 6 && count_c == 5 && count_n == 1 && count_o == 0 && count_s == 0 {
        return Some("pyridine".to_string());
    }
    // Recognize furan: 5-membered ring, 4 carbons and 1 oxygen.
    if ring.len() == 5 && count_c == 4 && count_o == 1 && count_n == 0 && count_s == 0 {
        return Some("furan".to_string());
    }
    // Recognize thiophene: 5-membered ring, 4 carbons and 1 sulfur.
    if ring.len() == 5 && count_c == 4 && count_s == 1 && count_o == 0 && count_n == 0 {
        return Some("thiophene".to_string());
    }
    // NEW: Recognize pyrrole or pyrrolidine:
    // For a 5-membered ring with 4 carbons and 1 nitrogen (and no oxygen or sulfur),
    // if the ring is aromatic we assume pyrrole; otherwise, we assume pyrrolidine.
    if ring.len() == 5 && count_c == 4 && count_n == 1 && count_o == 0 && count_s == 0 {
        if ring.iter().all(|&n| graph[n].is_aromatic()) {
            return Some("pyrrole".to_string());
        } else {
            return Some("pyrrolidine".to_string());
        }
    }
    // Fallback: build a generic name by concatenating element symbols and counts.
    let mut freq: HashMap<String, usize> = HashMap::new();
    for &n in ring {
        let sym = graph[n].symbol().to_string();
        *freq.entry(sym).or_insert(0) += 1;
    }
    let mut parts: Vec<String> = freq
        .into_iter()
        .map(|(s, c)| format!("{}{}", s, c))
        .collect();
    parts.sort();
    Some(parts.join(""))
}

/// Returns true if the given ring is imidazole (a five‐membered ring with three carbons and two nitrogens).
fn is_imidazole_ring(graph: &MoleculeGraph, ring: &[NodeIndex]) -> bool {
    if ring.len() != 5 {
        return false;
    }
    let mut count_c = 0;
    let mut count_n = 0;
    for &n in ring {
        let elem = graph[n];
        if elem == Element::C {
            count_c += 1;
        } else if elem == Element::N {
            count_n += 1;
        }
    }
    count_c == 3 && count_n == 2
}

/// Extend the heterocyclic substituent naming so that imidazole is handled specially.
fn heterocycle_substituent_name(graph: &MoleculeGraph, start: NodeIndex) -> Option<String> {
    trace!("Checking for heterocycle substituent at node {:?}", start);
    if let Some(ring) = find_ring_containing_node(graph, start) {
        trace!("Found ring containing substituent: {:?}", ring);
        // If the ring is imidazole, return the preferred substituent name.
        if is_imidazole_ring(graph, &ring) {
            // Note: Histidine is named as 2‑amino‑3‑(1H‑imidazol‑4‑yl)propanoic acid.
            // That is, the side chain substituent should be "1H-imidazol-4-yl"
            return Some("1H-imidazol-4-yl".to_string());
        }
        // Otherwise, fallback: use the ring's base name (as determined by base_name_for_ring)
        // and append a "yl" suffix.
        if let Some(base) = base_name_for_ring(graph, &ring) {
            // SPECIAL CASE: if the ring is benzene, we want "phenyl" not "benzyl".
            if base == "benzene" {
                return Some("phenyl".to_string());
            }
            if base.ends_with("ene") {
                let subst = base[..base.len()-3].to_string() + "yl";
                return Some(subst);
            }
            return Some(base + "yl");
        }
    }
    None
}


/// Adjusted naming for substituents:
/// - A double-bonded oxygen (to the parent) is named "oxo" (for aldehydes).
/// - A single-bonded oxygen that connects to another carbon is named as an ether (e.g. "methoxy").
/// - Otherwise, oxygen is "hydroxy".
fn name_substituent(graph: &MoleculeGraph, neighbor: NodeIndex, exclude: NodeIndex) -> String {
    // NEW: If the neighbor is part of a heterocycle, use our special handling.
    // info!("Name substituent: {:?}", neighbor);
    // info!("find_ring_containing_node {neighbor:?}: {:?}", find_ring_containing_node(graph, neighbor));
    // info!("find_ring_containing_node {neighbor:?}: {:?}", find_ring_containing_node(graph, neighbor));
    if let Some(ring) = find_ring_containing_node(graph, neighbor) {
        trace!("Found ring containing substituent: {:?}", ring);
        if !ring.iter().all(|&n| graph[n].is_carbon()) {
            if let Some(name) = heterocycle_substituent_name(graph, neighbor) {
                return name;
            }
        } else {
            // If the ring is fully carbon, fallback to the generic naming.
            if let Some(name) = base_name_for_ring(graph, &ring) {
                if name == "benzene" {
                    return "phenyl".to_string();
                }
                info!("Ring base name: {}", name);
                return name + "yl";
            }
        }
    } else {
        trace!("No ring found for substituent: {:?}", neighbor);
    }

    match graph[neighbor] {
        Element::C => {
            // NEW: Check if the carbon substituent fits a sulfhydryl pattern.
            if let Some(sh) = detect_sulfhydryl_substituent(graph, neighbor, exclude) {
                return sh;
            }
            // NEW: If the carbon substituent fits a carboxamide pattern, use "carbamoyl".
            if let Some(name) = detect_carboxamide_substituent(graph, neighbor, exclude) {
                trace!("Detected carboxamide substituent: {}", name);
                return name;
            } else {
                trace!("No carboxamide substituent detected");
            }
            if let Some(iso) = detect_iso_substituent(graph, neighbor, exclude) {
                return iso;
            }
            // For carbon branches, count how many connected carbon atoms are in this substituent.
            let branch_size = count_branch_carbons(graph, neighbor, exclude);
            // trace!("Chain substituent: {:?} in {graph:#?}", neighbor);

            match branch_size {
                1 => "methyl".to_string(),
                2 => "ethyl".to_string(),
                3 => "propyl".to_string(),
                4 => "butyl".to_string(),
                5 => "pentyl".to_string(),
                6 => "hexyl".to_string(),
                7 => "heptyl".to_string(),
                8 => "octyl".to_string(),
                9 => "nonyl".to_string(),
                10 => "decyl".to_string(),
                11 => "undecyl".to_string(),
                12 => "dodecyl".to_string(),
                _ => format!("{}-carbon alkyl", branch_size),
            }
        }
        // NEW: add an arm for sulfur.
        Element::S => {
            // We assume that a direct S attached as a substituent should be named "sulfhydryl".
            return "sulfhydryl".to_string();
        }
        Element::O => {
            // Check the bond type between the parent and oxygen.
            if let Some(edge) = graph.find_edge(exclude, neighbor) {
                if graph.edge_weight(edge) == Some(&Bond::Double) {
                    // A double bond indicates a carbonyl (aldehyde) group.
                    return "oxo".to_string();
                }
            }
            // For single-bonded oxygen, determine if it connects to another carbon.
            let carbon_neighbors: Vec<_> = graph
                .neighbors(neighbor)
                .filter(|&nbr| nbr != exclude && graph[nbr] == Element::C)
                .collect();
            if !carbon_neighbors.is_empty() {
                // If an oxygen is bridging to another carbon, treat it as an ether substituent.
                let branch_size = count_branch_carbons(graph, carbon_neighbors[0], neighbor);
                let alkoxy = match branch_size {
                    1 => "methoxy",
                    2 => "ethoxy",
                    3 => "propoxy",
                    4 => "butoxy",
                    5 => "pentoxy",
                    _ => "alkoxy",
                };
                return alkoxy.to_string();
            }
            // Default: if no additional carbon is found, treat as hydroxy.
            "hydroxy".to_string()
        }
        Element::N => {
            // NEW: Try to detect a guanidino group first.
            if let Some(guanidino) = detect_guanidino_substituent(graph, neighbor, exclude) {
                return guanidino;
            }
            if let Some(edge) = graph.find_edge(exclude, neighbor) {
                if let Some(bond) = graph.edge_weight(edge) {
                    if *bond == Bond::Triple {
                        return "cyano".to_string();
                    }
                }
            }
            "amino".to_string()
        }
        Element::F => "fluoro".to_string(),
        Element::Cl => "chloro".to_string(),
        Element::Br => "bromo".to_string(),
        Element::I => "iodo".to_string(),
        Element { kind, .. } if kind.is_r_group() => {
            // Handle R-groups as generic substituents.
            kind.symbol().to_string()
        }
        unknown => unknown.symbol().to_string(),
    }
}

fn count_branch_carbons(graph: &MoleculeGraph, start: NodeIndex, exclude: NodeIndex) -> usize {
    let mut visited = HashSet::new();
    visited.insert(exclude);
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if visited.insert(node) {
            for neighbor in graph.neighbors(node) {
                if !visited.contains(&neighbor) && graph[neighbor].is_carbon() {
                    stack.push(neighbor);
                }
            }
        }
    }
    visited.len() - 1
}

// =============================================================================
// Modified helper: use saturating_sub to avoid subtraction overflow
fn choose_numbering(main_chain: &[NodeIndex], substituents: &[(usize, String)]) -> NumberingOrder {
    if substituents.is_empty() {
        return NumberingOrder::Original;
    }
    let n = main_chain.len();
    let mut forward_locants: Vec<usize> = substituents.iter().map(|(pos, _)| *pos).collect();
    forward_locants.sort();
    let mut reverse_locants: Vec<usize> = substituents
        .iter()
        .map(|(pos, _)| (n + 1).saturating_sub(*pos))
        .collect();
    reverse_locants.sort();
    if reverse_locants < forward_locants {
        NumberingOrder::Reversed
    } else {
        NumberingOrder::Original
    }
}


fn assign_substituent_positions(
    final_main_chain: &Vec<NodeIndex>,
    substituents: &Vec<(usize, String)>,
    numbering_order: NumberingOrder,
) -> Vec<Substituent> {
    let n = final_main_chain.len();
    substituents
        .iter()
        .map(|(pos, name)| {
            let new_pos = match numbering_order {
                NumberingOrder::Original => *pos,
                NumberingOrder::Reversed => n + 1 - pos,
            };
            Substituent {
                position: new_pos,
                name: name.clone(),
            }
        })
        .collect()
}

fn format_substituents(substituents: &Vec<Substituent>) -> String {
    if substituents.is_empty() {
        return String::new();
    }
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for sub in substituents {
        groups.entry(sub.name.clone()).or_default().push(sub.position);
    }
    let mut parts = Vec::new();
    for (name, mut positions) in groups {
        positions.sort();
        // Join locants within the same group with commas (if there are multiple locants).
        let locant_str = positions
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let prefix = match positions.len() {
            1 => "",
            2 => "di",
            3 => "tri",
            4 => "tetra",
            _ => "",
        };
        // For each group, produce a string like "2-amino" or "5-guanidino"
        parts.push(format!("{}-{}{}", locant_str, prefix, name));
    }
    // Join the different substituent groups with hyphens.
    parts.join("-")
}

/// Returns the IUPAC name for a cyclic molecule (cycloalkane possibly with substituents).
fn iupac_cyclo_name(graph: &MoleculeGraph, ring: &Vec<NodeIndex>) -> String {
    // Check if the ring is fully aromatic.
    let is_aromatic_ring = ring
        .iter()
        .all(|&n| graph[n].is_aromatic() && graph[n].is_carbon());

    if is_aromatic_ring && ring.len() == 6 {
        // Identify substituents on the aromatic ring.
        let subs = identify_substituents_on_ring(graph, ring);
        if subs.is_empty() {
            return "benzene".to_string();
        } else {
            // Determine the best numbering for the ring by trying all rotations and directions.
            let (_best_ring, best_subs) = best_ring_numbering(ring, graph, &subs);
            let substituent_prefix = format_substituents(&best_subs);
            return format!("{}benzene", substituent_prefix);
        }
    }

    // Maybe its not aromatic, but has alternating double bonds.
    if ring.len() == 6 && ring.iter().all(|&n| graph[n].is_carbon()) {
        // Count the number of double bonds.
        let double_bonds = ring
            .windows(2)
            .filter(|&pair| {
                if let Some(edge) = graph.find_edge(pair[0], pair[1]) {
                    graph.edge_weight(edge) == Some(&Bond::Double)
                } else {
                    false
                }
            })
            .count();

        // Confirm no double bonds are adjacent.
        if double_bonds == 3 {
            // Check for alternating double bonds.
            let alternating = ring
                .windows(2)
                .enumerate()
                .all(|(i, pair)| {
                    if let Some(edge) = graph.find_edge(pair[0], pair[1]) {
                        if i % 2 == 0 {
                            graph.edge_weight(edge) == Some(&Bond::Double)
                        } else {
                            graph.edge_weight(edge) == Some(&Bond::Single)
                        }
                    } else {
                        false
                    }
                });
            if alternating {
                // Identify substituents on the aromatic ring.
                let subs = identify_substituents_on_ring(graph, ring);
                if subs.is_empty() {
                    return "benzene".to_string();
                } else {
                    // Determine the best numbering for the ring by trying all rotations and directions.
                    let (_best_ring, best_subs) = best_ring_numbering(ring, graph, &subs);
                    let substituent_prefix = format_substituents(&best_subs);
                    return format!("{}benzene", substituent_prefix);
                }
            }
        }
    }

    // Non-aromatic (or aromatic but not 6-membered) cyclic naming.
    let subs = identify_substituents_on_ring(graph, ring);
    info!("Substituents on ring: {:?}", subs);
    let (best_ring, best_subs) = best_ring_numbering(ring, graph, &subs);
    let substituent_prefix = format_substituents(&best_subs);

    // Determine the base name for the ring.
    let n = best_ring.len();
    let base_name = match n {
        3 => "cyclopropane",
        4 => "cyclobutane",
        5 => "cyclopentane",
        6 => "cyclohexane",
        7 => "cycloheptane",
        8 => "cyclooctane",
        9 => "cyclononane",
        10 => "cyclodecane",
        _ => return format!("{}-membered ring (not supported)", n),
    };

    if substituent_prefix.is_empty() {
        base_name.to_string()
    } else {
        format!("{}{}", substituent_prefix, base_name)
    }
}

/// For a cyclic structure the numbering is ambiguous due to rotations and directions.
/// Try all rotations (and reversed order) to choose the numbering that minimizes
/// the sorted list of substituent locants.
fn best_ring_numbering(
    ring: &Vec<NodeIndex>,
    graph: &MoleculeGraph,
    subs: &Vec<(usize, String)>,
) -> (Vec<NodeIndex>, Vec<Substituent>) {
    let n = ring.len();
    let mut best_order = ring.clone();
    let mut best_subs: Vec<Substituent> = subs
        .iter()
        .map(|(pos, name)| Substituent {
            position: *pos,
            name: name.clone(),
        })
        .collect();
    best_subs.sort_by_key(|s| s.position);
    let mut best_locants: Vec<usize> = best_subs.iter().map(|s| s.position).collect();

    // Generate all rotations and directions.
    for start in 0..n {
        let candidate = ring[start..]
            .iter()
            .chain(ring[..start].iter())
            .cloned()
            .collect::<Vec<_>>();
        for &reverse in &[false, true] {
            let candidate_order = if reverse {
                let mut rev = candidate.clone();
                rev.reverse();
                rev
            } else {
                candidate.clone()
            };
            // Reassign substituents positions based on candidate_order.
            let candidate_subs = assign_substituents_on_ring(&candidate_order, graph);
            let mut candidate_locants: Vec<usize> =
                candidate_subs.iter().map(|s| s.position).collect();
            candidate_locants.sort();
            if candidate_locants < best_locants {
                best_locants = candidate_locants;
                best_order = candidate_order;
                best_subs = candidate_subs;
            }
        }
    }
    (best_order, best_subs)
}

/// Given a ring ordering (a vector of NodeIndex in cyclic order),
/// assign substituent positions based on that order.
fn assign_substituents_on_ring(
    ring_order: &Vec<NodeIndex>,
    graph: &MoleculeGraph,
) -> Vec<Substituent> {
    let mut subs = Vec::new();
    let ring_set: HashSet<_> = ring_order.iter().cloned().collect();
    for (i, &ring_node) in ring_order.iter().enumerate() {
        for neighbor in graph.neighbors(ring_node) {
            if !ring_set.contains(&neighbor) {
                let sub_name = name_substituent(graph, neighbor, ring_node);
                subs.push(Substituent {
                    position: i + 1, // positions are 1-indexed
                    name: sub_name,
                });
            }
        }
    }
    subs
}

fn is_in_ring(graph: &MoleculeGraph, node: NodeIndex) -> bool {
    find_ring_containing_node(graph, node).is_some()
}
// --- NEW: Acyl Chain Determination ---
// Starting at the acyl carbon, find the longest path that follows only carbon atoms.
// This chain represents the acid-derived portion of the ester.
fn find_acyl_chain(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> Vec<NodeIndex> {
    let mut visited = HashSet::new();
    dfs_longest_path(graph, acyl_carbon, &mut visited)
}

fn dfs_longest_nonring_path(graph: &MoleculeGraph, current: NodeIndex, visited: &mut HashSet<NodeIndex>) -> Vec<NodeIndex> {
    visited.insert(current);
    let mut longest = vec![current];
    for nbr in graph.neighbors(current) {
        if graph[nbr] == Element::C && !is_in_ring(graph, nbr) && !visited.contains(&nbr) {
            let path = dfs_longest_nonring_path(graph, nbr, visited);
            if 1 + path.len() > longest.len() {
                let mut candidate = vec![current];
                candidate.extend(path);
                longest = candidate;
            }
        }
    }
    longest
}

// --- Existing DFS for longest carbon path ---
fn dfs_longest_path(
    graph: &MoleculeGraph,
    current: NodeIndex,
    visited: &mut HashSet<NodeIndex>,
) -> Vec<NodeIndex> {
    visited.insert(current);
    let mut longest = vec![current];
    for neighbor in graph.neighbors(current) {
        if graph[neighbor] == Element::C && !visited.contains(&neighbor) {
            let mut new_visited = visited.clone();
            let path = dfs_longest_path(graph, neighbor, &mut new_visited);
            if 1 + path.len() > longest.len() {
                let mut candidate = vec![current];
                candidate.extend(path);
                longest = candidate;
            }
        }
    }
    longest
}

/// For an alkoxy substituent, if it is part of an aromatic ring, return its substituent name.
/// For example, benzene becomes "phenyl". In general, if the ring base name ends with "ene",
/// replace that ending with "yl". Otherwise, simply add "yl".
fn alkoxy_ring_name(graph: &MoleculeGraph, start: NodeIndex) -> Option<String> {
    if let Some(sub_ring) = find_ring_containing_node(graph, start) {
        // Check that the ring is fully aromatic.
        if sub_ring.iter().all(|&n| graph[n].is_carbon()) {
            if let Some(base) = base_name_for_ring(graph, &sub_ring) {
                if base == "benzene" {
                    return Some("phenyl".to_string());
                }
                if base.ends_with("ene") {
                    let subst = base[..base.len() - 3].to_string() + "yl";
                    return Some(subst);
                }
                return Some(base + "yl");
            }
        }
    }
    None
}

/// Returns the conventional locant for the attachment in a heterocycle.
/// This function computes both the clockwise and anticlockwise distances from
/// the first heteroatom (the first non‑carbon in best_order) and returns the smaller distance plus one.
fn conventional_locant_for_heterocycle(
    graph: &MoleculeGraph,
    best_order: &[NodeIndex],
    attach: NodeIndex,
) -> usize {
    let n = best_order.len();
    let h = best_order
        .iter()
        .position(|&node| graph[node].is_heteroatom())
        .unwrap_or(0);
    let a = best_order
        .iter()
        .position(|&node| node == attach)
        .unwrap_or(0);
    let d1 = if a >= h { a - h } else { a + n - h };
    let d2 = n - d1;
    d1.min(d2) + 1
}
/// New version of the acyl‐group naming function for ester contexts.
/// This function “knows” it must produce an ester (–oate or –carboxylate) name.
fn name_acyl_group_ester(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> String {
    // If the acyl carbon is directly bonded to an aromatic carbon in a six‐membered all‐carbon ring,
    // assume benzoate.
    let is_benzo = graph.neighbors(acyl_carbon).any(|nbr| {
        if graph[nbr].is_carbon() {
            if let Some(ring) = find_ring(graph) {
                ring.contains(&nbr) && ring.len() == 6 && ring.iter().all(|&n| graph[n].is_carbon())
            } else {
                false
            }
        } else {
            false
        }
    });
    if is_benzo {
        return "benzoate".to_string();
    }

    // If the acyl carbon is part of a ring (i.e. the acid fragment is heterocyclic),
    // use our ester‐aware heterocyclic naming:
    if let Some(ring) = find_ring_containing_node(graph, acyl_carbon) {
        return name_heterocyclic_acid_ester(graph, &ring, acyl_carbon);
    }

    // If the acyl carbon is exocyclic to a heterocycle,
    // try our ester‐aware exocyclic naming:
    if let Some(exo_name) = name_exocyclic_acid_ester(graph, acyl_carbon) {
        return exo_name;
    }

    // Otherwise, fallback to the acyclic naming branch.
    let acyl_chain = find_acyl_chain(graph, acyl_carbon);
    // (Keep your special-case for oxalic acid, etc.)
    // Also, incorporate substituent naming as before.
    // Finally, form the name using an “oate” suffix:
    let subs = identify_substituents(graph, &acyl_chain)
        .into_iter().filter(|(pos, _)| *pos != 1)
        .collect::<Vec<_>>();
    let numbering_order = choose_numbering(&acyl_chain, &subs);
    let final_chain = match numbering_order {
        NumberingOrder::Original => acyl_chain.clone(),
        NumberingOrder::Reversed => {
            let mut rev = acyl_chain.clone();
            rev.reverse();
            rev
        }
    };
    let subs_final = assign_substituent_positions(&final_chain, &subs, numbering_order);
    let substituent_prefix = format_substituents(&subs_final);
    let chain_length = final_chain.len();
    let root = match chain_length {
        1 => "methan",
        2 => "ethan",
        3 => "propan",
        4 => "butan",
        5 => "pentan",
        6 => "hexan",
        7 => "heptan",
        8 => "octan",
        9 => "nonan",
        10 => "decan",
        n => return format!("{}-carbon acid chain (not fully supported)", n),
    };
    format!("{}oate", format!("{}{}", substituent_prefix.trim_end_matches('-'), root))
}

/// Ester‐aware exocyclic acid naming.
/// Instead of returning a “‑carboxaldehyde” or “‑carboxylic acid” name,
/// it returns a name ending in “‑carboxylate.”
fn name_exocyclic_acid_ester(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> Option<String> {
    for nbr in graph.neighbors(acyl_carbon) {
        if graph[nbr] == Element::O { continue; }
        if let Some(ring) = find_ring_containing_node(graph, nbr) {
            if let Some(base) = base_name_for_ring(graph, &ring) {
                if base != "benzene" {
                    let (best_order, _) = best_heterocycle_numbering(&ring, graph);
                    let locant = conventional_locant_for_heterocycle(graph, &best_order, nbr);
                    return Some(format!("{}-{}-carboxylate", base, locant));
                }
            }
        }
    }
    None
}


/// Ester‐aware heterocyclic acid naming.
/// This version is used when the acyl carbon is part of (or directly attached to) a heterocycle
/// and should return the “–carboxylate” suffix.
fn name_heterocyclic_acid_ester(graph: &MoleculeGraph, ring: &[NodeIndex], acyl_carbon: NodeIndex) -> String {
    if let Some(base) = base_name_for_ring(graph, ring) {
        let (best_order, _) = best_heterocycle_numbering(ring, graph);
        if let Some(pos) = best_order.iter().position(|&n| n == acyl_carbon) {
            let locant = pos + 1;
            if base == "benzene" {
                return "benzenecarboxylate".to_string();
            }
            return format!("{}-{}-carboxylate", base, locant);
        }
        return format!("{}carboxylate", base);
    }
    "carboxylate".to_string()
}

fn iupac_single_ester_name(
    graph: &MoleculeGraph,
    acyl_carbon: NodeIndex,
    alkoxy_oxygen: NodeIndex,
) -> String {
    debug!("Naming single ester: acyl carbon {:?}, alkoxy oxygen {:?}", acyl_carbon, alkoxy_oxygen);
    // Use the ester‐aware acyl group naming.
    let acyl_name = name_acyl_group_ester(graph, acyl_carbon);
    if let Some(alkoxy_carb) = graph
        .neighbors(alkoxy_oxygen)
        .find(|&nbr| nbr != acyl_carbon && graph[nbr].is_carbon())
    {
        if let Some(ring_name) = alkoxy_ring_name(graph, alkoxy_carb) {
            return format!("{}-{}", ring_name, acyl_name);
        }
        if is_tert_butyl(graph, alkoxy_carb, alkoxy_oxygen) {
            return format!("tert-butyl-{}", acyl_name);
        }
        if let Some(iso) = detect_iso_substituent(graph, alkoxy_carb, alkoxy_oxygen) {
            return format!("{}-{}", iso, acyl_name);
        }
        let branch_size = count_branch_carbons(graph, alkoxy_carb, alkoxy_oxygen);
        let alkoxy_name = match branch_size {
            1 => "methyl",
            2 => "ethyl",
            3 => "propyl",
            4 => "butyl",
            5 => "pentyl",
            _ => "alkyl",
        };
        return format!("{}-{}", alkoxy_name, acyl_name);
    }
    acyl_name
}

/// Given an acyl (carbonyl) carbon that is not in a ring, check its non‑oxygen neighbors.
/// If one is in a heterocyclic ring (i.e. the acid fragment is exocyclic),
/// then use conventional numbering to form a name like "<base>-<locant>-carboxylate".
fn name_exocyclic_acid(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> Option<String> {
    for nbr in graph.neighbors(acyl_carbon) {
        if graph[nbr] == Element::O {
            continue;
        }
        if let Some(ring) = find_ring_containing_node(graph, nbr) {
            if let Some(base) = base_name_for_ring(graph, &ring) {
                if base != "benzene" {
                    // only treat heterocycles specially
                    let (best_order, _) = best_heterocycle_numbering(&ring, graph);
                    let locant = conventional_locant_for_heterocycle(graph, &best_order, nbr);
                    return Some(format!("{}-{}-carboxylate", base, locant));
                }
            }
        }
    }
    None
}

/// If the aldehyde (carbonyl) group is attached exocyclically to a heterocycle,
/// this function computes the conventional locant from the heterocycle’s best ordering
/// and returns a name like "furan-2-carboxaldehyde".
fn name_exocyclic_aldehyde(graph: &MoleculeGraph, carbonyl_carbon: NodeIndex) -> Option<String> {
    for nbr in graph.neighbors(carbonyl_carbon) {
        if graph[nbr] == Element::O {
            continue;
        }
        if let Some(ring) = find_ring_containing_node(graph, nbr) {
            if let Some(base) = base_name_for_ring(graph, &ring) {
                if base != "benzene" {
                    // only for heterocycles
                    let (best_order, _) = best_heterocycle_numbering(&ring, graph);
                    let locant = conventional_locant_for_heterocycle(graph, &best_order, nbr);
                    return Some(format!("{}-{}-carboxaldehyde", base, locant));
                }
            }
        }
    }
    None
}

/// Ester‐or free–acid heterocyclic naming helper for amino acids.
/// If the acid carbon is in the ring, use its locant; if it’s exocyclic, use the locant
/// of the ring atom it attaches to.
fn name_heterocyclic_acid_free(
    graph: &MoleculeGraph,
    ring: &[NodeIndex],
    acyl_carbon: NodeIndex,
) -> String {
    if let Some(base) = base_name_for_ring(graph, ring) {
        trace!("Base name for ring: {}", base);
        // Determine the best numbering for the ring.
        let (best_order, best_hetero) = best_heterocycle_numbering(ring, graph);
        trace!("Best heterocycle numbering: {:?}", best_order);
        trace!("Best heteroatom: {:?}", best_hetero);
        
        // First try: if the acid carbon is already in the ring, use its position.
        if let Some(pos) = best_order.iter().position(|&n| n == acyl_carbon) {
            let locant = pos + 1;
            return format!("{}-{}-carboxylic acid", base, locant);
        }
        // Otherwise, the acid carbon is exocyclic.
        // Find which atom in the ring is attached to it and use that atom's position.
        for nbr in graph.neighbors(acyl_carbon) {
            if ring.contains(&nbr) {
                if let Some(pos) = best_order.iter().position(|&n| n == nbr) {
                    let locant = pos + 1;
                    return format!("{}-{}-carboxylic acid", base, locant);
                }
            }
        }
        // If no attached ring atom is found (should not happen for amino acids), fall back.
        return format!("{}-carboxylic acid", base);
    }
    "carboxylic acid".to_string()
}

/// Given a ring that contains the acyl carbon, returns a systematic acid name
/// of the form "<base>-<locant>-carboxylate". The numbering is determined by the
/// existing best_heterocycle_numbering routine.
fn name_heterocyclic_acid(
    graph: &MoleculeGraph,
    ring: &[NodeIndex],
    acyl_carbon: NodeIndex,
) -> String {
    if let Some(base) = base_name_for_ring(graph, ring) {
        let (best_order, _) = best_heterocycle_numbering(ring, graph);
        if let Some(pos) = best_order.iter().position(|&n| n == acyl_carbon) {
            let locant = pos + 1;
            // For benzene, return a conventional name.
            if base == "benzene" {
                return "benzenecarboxylate".to_string();
            }
            return format!("{}-{}-carboxylate", base, locant);
        }
        return format!("{}carboxylate", base);
    }
    "carboxylate".to_string()
}

/// Helper: returns a ring (if any) that contains the given node.
fn find_ring_containing_node(graph: &MoleculeGraph, node: NodeIndex) -> Option<Vec<NodeIndex>> {
    // if let Some(ring) = find_ring(graph) {
    //     if ring.contains(&node) {
    //         return Some(ring);
    //     }
    // }
    // None
    
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    info!("Checking for ring containing node: {node:?}");
    for nbr in graph.neighbors(node) {
        visited.clear();
        stack.clear();
        // visited.insert(node);
        // stack.push(node);
        info!("Starting DFS from neighbor: {nbr:?}");
        if let Some(mut cycle) = dfs_find_cycle_generic(graph, Some(node), nbr, None, &mut stack, &mut visited) {
            if cycle.contains(&node) {
                cycle.sort_by_key(|&n| n.index());
                return Some(cycle);
            }
        }
    }

    // warn!("No cycle found containing {node:?}");
    None

    // let mut visited = HashSet::new();
    // let mut stack = Vec::new();
    // if let Some(cycle) = dfs_find_cycle_generic(graph, node, None, &mut stack, &mut visited) {
    //     debug!("Found cycle containing {node:?}: {:?}", cycle);
    //     return Some(cycle);
    // }
    // warn!("No cycle found containing {node:?}");
    // None

    // Perform a DFS to find a cycle containing the node.
    // let mut visited = HashSet::new();
    // let mut stack = Vec::new();
    // if graph[node].is_carbon() {
    //     for nbr in graph.neighbors(node) {
    //         if graph[nbr].is_carbon() && !visited.contains(&nbr) {
    //             if let Some(cycle) = dfs_find_cycle(graph, nbr, None, &mut stack, &mut visited) {
    //                 if cycle.contains(&node) {
    //                     return Some(cycle);
    //                 }
    //             }
    //         }
    //         visited.clear();
    //         stack.clear();
    //     }
    // }
}

/// --- Alkoxy Naming Helpers ---
/// Detects whether the candidate alkoxy carbon is a tert-butyl group.
fn is_tert_butyl(graph: &MoleculeGraph, node: NodeIndex, exclude: NodeIndex) -> bool {
    if !graph[node].is_carbon() {
        return false;
    }
    let carbon_neighbors: Vec<NodeIndex> = graph
        .neighbors(node)
        .filter(|&nbr| nbr != exclude && graph[nbr].is_carbon())
        .collect();
    if carbon_neighbors.len() != 3 {
        return false;
    }
    // Check that each neighbor is terminal.
    for nbr in carbon_neighbors {
        let count = graph
            .neighbors(nbr)
            .filter(|&x| x != node && graph[x].is_carbon())
            .count();
        if count != 0 {
            return false;
        }
    }
    true
}

// Returns an “iso” name for a substituent if it fits a recognized pattern.
/// Currently supports detection of isopropyl (3‑carbon) and isobutyl (4‑carbon) groups.
fn detect_iso_substituent(
    graph: &MoleculeGraph,
    start: NodeIndex,
    exclude: NodeIndex,
) -> Option<String> {
    // For isopropyl:
    // The attached carbon (start) should have exactly 2 carbon neighbors (besides exclude),
    // and each such neighbor should be terminal (i.e. no further carbon neighbors besides start).
    let carb_neighbors: Vec<_> = graph
        .neighbors(start)
        .filter(|&n| n != exclude && graph[n].is_carbon())
        .collect();
    if carb_neighbors.len() == 2 {
        if carb_neighbors.iter().all(|&nbr| {
            graph
                .neighbors(nbr)
                .filter(|&n| n != start && graph[n].is_carbon())
                .count()
                == 0
        }) {
            return Some("isopropyl".to_string());
        }
    }
    // For isobutyl:
    // The attached carbon (start) should be primary (exactly one carbon neighbor besides exclude),
    // and its unique neighbor should have exactly 2 terminal carbon neighbors.
    let carb_neighbors: Vec<_> = graph
        .neighbors(start)
        .filter(|&n| n != exclude && graph[n].is_carbon())
        .collect();
    if carb_neighbors.len() == 1 {
        let second = carb_neighbors[0];
        let second_neighbors: Vec<_> = graph
            .neighbors(second)
            .filter(|&n| n != start && graph[n].is_carbon())
            .collect();
        if second_neighbors.len() == 2
            && second_neighbors.iter().all(|&nbr| {
                graph
                    .neighbors(nbr)
                    .filter(|&n| n != second && graph[n].is_carbon())
                    .count()
                    == 0
            })
        {
            return Some("isobutyl".to_string());
        }
    }
    None
}

fn detect_acid_global(graph: &MoleculeGraph) -> Option<(Vec<NodeIndex>, usize)> {
    debug!("Running detect_acid_global");
    let mut candidate_chains = Vec::new();
    // Iterate over all carbons looking for a carboxyl carbon.
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            let mut has_dbl = false;
            let mut has_sng = false;
            for nbr in graph.neighbors(node) {
                if graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        match graph.edge_weight(edge) {
                            Some(Bond::Double) => has_dbl = true,
                            Some(Bond::Single) => has_sng = true,
                            _ => {}
                        }
                    }
                }
            }
            if has_dbl && has_sng {
                // This node qualifies as a carboxyl carbon.
                // Find the longest carbon chain (using only C–C bonds) starting at this node.
                let chain = dfs_longest_path(graph, node, &mut HashSet::new());
                candidate_chains.push(chain);
            }
        }
    }
    if candidate_chains.is_empty() {
        debug!("No acid carbon found");
        return None;
    }
    // Choose the longest chain among the candidates.
    let longest = candidate_chains.into_iter().max_by_key(|chain| chain.len())?;
    // Return this chain with acid index 0 (i.e. the acid carbon is first).
    Some((longest, 0))
}

/// Detect a carboxamide substituent pattern.
/// If the substituent branch (starting at the carbon attached to the main chain)
/// has a double‐bonded oxygen and a single‐bonded nitrogen (i.e. fits –C(N)=O),
/// then we interpret it as a carboxamide group and return "carbamoyl".
fn detect_carboxamide_substituent(
    graph: &MoleculeGraph,
    start: NodeIndex,
    exclude: NodeIndex,
) -> Option<String> {
    trace!("Running detect_carboxamide_substituent");
    // Instead of looking at a neighbor-of-a-neighbor, check the "start" node itself.
    let mut has_oxo = false;
    let mut has_amino = false;
    for nbr in graph.neighbors(start) {
        if nbr == exclude {
            continue;
        }
        match graph[nbr] {
            Element::O => {
                if let Some(edge) = graph.find_edge(start, nbr) {
                    if graph.edge_weight(edge) == Some(&Bond::Double) {
                        has_oxo = true;
                    }
                }
            }
            Element::N => {
                if let Some(edge) = graph.find_edge(start, nbr) {
                    if graph.edge_weight(edge) == Some(&Bond::Single) {
                        has_amino = true;
                    }
                }
            }
            _ => {}
        }
    }
    if has_oxo && has_amino {
        return Some("carbamoyl".to_string());
    }
    trace!("No carboxamide substituent detected");
    None
}

fn detect_amino_acid(graph: &MoleculeGraph) -> Option<String> {
    trace!("Running detect_amino_acid");

    // Locate the acid (carboxyl) carbon.
    let acid = find_acid_carbon(graph)?;
    // First, if the acid carbon is directly in a ring, use the heterocyclic free–acid naming.
    // if let Some(ring) = find_ring_containing_node(graph, acid) {
    //     trace!("Acid carbon is in a ring: {:?}", ring);
    //     return Some(name_heterocyclic_acid_free(graph, &ring, acid));
    // }
    // Otherwise, check if any non‐oxygen neighbor of the acid carbon belongs to a ring.
    for nbr in graph.neighbors(acid) {
        if graph[nbr] != Element::O {
            if let Some(ring) = find_ring_containing_node(graph, nbr) {
                trace!("Acid carbon is exocyclic to a ring: {:?}", ring);
                // Here we “borrow” the ring from the neighbor and treat the acid carbon as exocyclic.
                return Some(name_heterocyclic_acid_free(graph, &ring, acid));
            }
        }
    }
    // If no attached ring is found, fallback to the linear (acyclic) amino–acid naming.
    let chain = find_acid_chain(graph, acid);
    if chain.len() < 3 {
        return None;
    }

    // We want a 3-carbon backbone for an amino acid.
    // If the chain is exactly 3 carbons, we use it.
    // If it is 4 carbons long and the extra carbon (chain[3]) fits a carboxamide pattern
    // attached to the β–carbon (chain[2]), then we trim the backbone to 3 carbons.
    let backbone: Vec<NodeIndex>;
    let mut extra_subs: Vec<(usize, String)> = Vec::new();
    if chain.len() == 3 {
        backbone = chain.clone();
    } else if chain.len() == 4 {
        if let Some(carbamoyl) = detect_carboxamide_substituent(graph, chain[3], chain[2]) {
            // Use only the first three carbons as the backbone.
            backbone = vec![chain[0], chain[1], chain[2]];
            // Record the overriding substituent "carbamoyl" on carbon 3.
            extra_subs.push((3, carbamoyl));
        } else {
            return None;
        }
    } else {
        return None;
    }
    // Check that the α–carbon (backbone[1]) carries an amino substituent.
    let alpha_subs = identify_substituents(graph, &[backbone[1]]);
    let mut has_amino = false;
    for (_, name) in alpha_subs.iter() {
        if name.contains("amino") {
            has_amino = true;
            break;
        }
    }
    if !has_amino {
        return None;
    }
    // Now, collect substituents on positions 2 and 3 from the backbone.
    let mut subs = identify_substituents(graph, &backbone);
    // Always ignore substituents on the acid carbon (position 1).
    subs.retain(|(pos, _)| *pos != 1);
    // If the backbone was derived from a 4–carbon chain (and trimmed),
    // then remove substituents on position 3 so that only the extra substituent remains.
    if chain.len() == 4 {
        subs.retain(|(pos, _)| *pos != 3);
    }
    // Then add the extra carbamoyl substituent.
    subs.extend(extra_subs);
    let prefix = format_substituents(&assign_substituent_positions(
        &backbone,
        &subs,
        NumberingOrder::Original,
    ));
    let root = "propan";
    Some(format!("{}oic acid", format!("{}{}", prefix.trim_end_matches('-'), root)))
}

/// Locate a carbon that qualifies as the carboxyl (acid) carbon.
/// It must be bonded to one oxygen with a double bond (carbonyl) and one oxygen with a single bond (–OH).
fn find_acid_carbon(graph: &MoleculeGraph) -> Option<NodeIndex> {
    for node in graph.node_indices() {
        if graph[node] == Element::C {
            let mut has_dbl = false;
            let mut has_sng = false;
            for nbr in graph.neighbors(node) {
                if graph[nbr] == Element::O {
                    if let Some(edge) = graph.find_edge(node, nbr) {
                        match graph.edge_weight(edge) {
                            Some(Bond::Double) => has_dbl = true,
                            Some(Bond::Single) => has_sng = true,
                            _ => {}
                        }
                    }
                }
            }
            if has_dbl && has_sng {
                return Some(node);
            }
        }
    }
    None
}

/// Starting from the given acid carbon, perform a DFS (following only C–C bonds)
/// to produce a chain. (For simplicity we reuse our DFS for the longest chain.)
fn find_acid_chain(graph: &MoleculeGraph, acid: NodeIndex) -> Vec<NodeIndex> {
    dfs_longest_path(graph, acid, &mut HashSet::new())
}
/// Returns true if the given node qualifies as an acid carbon 
/// (bonded to at least one double‐bonded O and one single‐bonded O).
fn is_acid_carbon(graph: &MoleculeGraph, node: NodeIndex) -> bool {
    let mut has_dbl = false;
    let mut has_sng = false;
    for nbr in graph.neighbors(node) {
        if graph[nbr] == Element::O {
            if let Some(edge) = graph.find_edge(node, nbr) {
                match graph.edge_weight(edge) {
                    Some(Bond::Double) => has_dbl = true,
                    Some(Bond::Single) => has_sng = true,
                    _ => {}
                }
            }
        }
    }
    has_dbl && has_sng
}

/// In the diacid case we want to ensure that the alkane root (e.g. "butan")
/// becomes "butane" (if it does not already end with an "e").
fn adjust_diacid_root(root: &str) -> String {
    if root.ends_with('e') {
        root.to_string()
    } else {
        format!("{}e", root)
    }
}

/// A modified free–acid naming routine that handles monoacid,
/// diacid and polyacid cases.
fn iupac_acid_name(graph: &MoleculeGraph) -> String {
    // First, locate the acid carbon.
    let acid_carbon = match find_acid_carbon(graph) {
        Some(c) => c,
        None => return "Unknown acid".to_string(),
    };
    // If the acid carbon is part of a ring, use the free–acid heterocyclic naming:
    if let Some(ring) = find_ring_containing_node(graph, acid_carbon) {
        trace!("Acid carbon is in a ring: {:?}", ring);
        return name_heterocyclic_acid_free(graph, &ring, acid_carbon);
    }
    // Build the acid chain (with acid carbon at the beginning).
    let mut chain = find_acid_chain(graph, acid_carbon);
    if chain[0] != acid_carbon {
        if let Some(pos) = chain.iter().position(|&n| n == acid_carbon) {
            chain.rotate_left(pos);
        }
    }
    // Count how many carbons in the chain are acid carbons.
    let acid_count = chain.iter().filter(|&&n| is_acid_carbon(graph, n)).count();
    // Identify substituents on the chain.
    let mut subs = identify_substituents(graph, &chain);
    // Remove substituents on the first carbon.
    subs.retain(|(pos, _)| *pos != 1);
    // If the chain is a diacid (both first and last carbons are acid carbons),
    // also remove substituents on the last carbon.
    let is_diacid = is_acid_carbon(graph, chain[0]) && is_acid_carbon(graph, *chain.last().unwrap());
    if is_diacid {
        subs.retain(|(pos, _)| *pos != chain.len());
    }
    // Force original numbering for acids.
    let numbering_order = NumberingOrder::Original;
    let final_chain = chain.clone();
    let subs_final = assign_substituent_positions(&final_chain, &subs, numbering_order);
    let substituent_prefix = format_substituents(&subs_final);
    
    // Determine the alkane root from chain length.
    let chain_length = final_chain.len();
    let root = match chain_length {
        1 => "methan",
        2 => "ethan",
        3 => "propan",
        4 => "butan",
        5 => "pentan",
        6 => "hexan",
        7 => "heptan",
        8 => "octan",
        9 => "nonan",
        10 => "decan",
        n => return format!("{}-carbon acid chain (not fully supported)", n),
    };

    // Now choose the final name based on the number of acid carbons.
    if acid_count == 2 {
        let diacid_root = adjust_diacid_root(root);
        format!("{}dioic acid", format!("{}{}", substituent_prefix.trim_end_matches('-'), diacid_root))
    } else if acid_count > 2 {
        let multiplier = match acid_count {
            3 => "tri",
            4 => "tetra",
            5 => "penta",
            6 => "hexa",
            7 => "hepta",
            8 => "octa",
            9 => "nona",
            10 => "deca",
            _ => "",
        };
        format!("{}{}-{}carboxylic acid", substituent_prefix.trim_end_matches('-'), root, multiplier)
    } else {
        format!("{}oic acid", format!("{}{}", substituent_prefix.trim_end_matches('-'), root))
    }
}

/// Returns the longest chain of carbons that are not part of any ring.
/// If none is found, returns an empty vector.
fn find_longest_nonring_carbon_chain(graph: &MoleculeGraph) -> Vec<NodeIndex> {
    let mut longest = Vec::new();
    for node in graph.node_indices() {
        if graph[node] == Element::C && !is_in_ring(graph, node) {
            let mut visited = HashSet::new();
            let chain = dfs_longest_nonring_path(graph, node, &mut visited);
            if chain.len() > longest.len() {
                longest = chain;
            }
        }
    }
    longest
}

/// Determines whether the given ring is attached to the acyclic parent chain
/// at only a single atom (and thus should be treated as a substituent rather than
/// the main parent structure).
fn is_ring_substituent(
    graph: &MoleculeGraph,
    ring: &Vec<NodeIndex>,
    acyclic_chain: &Vec<NodeIndex>,
) -> bool {
    let chain_set: HashSet<NodeIndex> = acyclic_chain.iter().cloned().collect();
    let common: Vec<NodeIndex> = ring.iter().filter(|n| chain_set.contains(n)).cloned().collect();
    // If the ring touches the chain at exactly one atom, then it is a substituent.
    common.len() == 1
}

/// Now, in the main naming routine we let the amino–acid branch override only if it applies,
/// and otherwise the generic acid branch (which now distinguishes diacids) takes over.
pub fn iupac_name(graph: &MoleculeGraph) -> String {
    debug!("Starting IUPAC naming process");
    if detect_ester_global(graph).is_some() {
        debug!("Ester detected globally; using acyclic naming branch");
        return iupac_acyclic_name(graph);
    }
    // Prioritize amino–acid naming if the pattern is detected.
    if let Some(name) = detect_amino_acid(graph) {
        debug!("Amino acid pattern detected; using amino acid naming branch");
        return name;
    } else {
        debug!("No amino acid pattern detected");
    }
    // Prioritize exocyclic aldehyde detection if present.
    if let Some(name) = detect_exocyclic_aldehyde(graph) {
        debug!("Exocyclic aldehyde detected; using that naming branch");
        return name;
    } else {
        debug!("No exocyclic aldehyde detected");
    }


    // Gather candidate parent structures.
    let acyclic_chain = find_longest_nonring_carbon_chain(graph);
    let ring_opt = find_ring(graph);

    // Decide the parent based on which has the higher priority.
    if !acyclic_chain.is_empty() {
        // If a nonring chain exists, consider the possibility that the ring is just a substituent.
        if let Some(ring) = ring_opt {
            // Rule-of-thumb: if the acyclic chain is longer than the ring,
            // then use the chain as the parent.
            if acyclic_chain.len() > ring.len() {
                debug!("Acyclic chain is longer than ring; using acyclic naming");
                return iupac_acyclic_name(graph);
            } else {
                // If lengths are comparable, check if the ring is attached in a manner
                // that suggests it is a substituent. For example, if the ring touches the chain at only one atom.
                if is_ring_substituent(graph, &ring, &acyclic_chain) {
                    debug!("Ring appears to be a substituent; using acyclic naming (phenyl groups will be named correctly)");
                    return iupac_acyclic_name(graph);
                } else {
                    debug!("Ring is the major component; using ring naming");
                    if ring.iter().all(|&n| graph[n].is_carbon()) {
                        return iupac_cyclo_name(graph, &ring);
                    } else {
                        return iupac_heterocyclo_name(graph, &ring);
                    }
                }
            }
        } else {
            debug!("No ring detected; using acyclic naming");
            return iupac_acyclic_name(graph);
        }
    } else if let Some(ring) = ring_opt {
        debug!("No nonring chain found; using ring naming");
        if ring.iter().all(|&n| graph[n].is_carbon()) {
            return iupac_cyclo_name(graph, &ring);
        } else {
            return iupac_heterocyclo_name(graph, &ring);
        }
    }

    debug!("Fallback to acyclic naming");
    iupac_acyclic_name(graph)
}

/// --- Existing DFS and Naming Routines (unchanged) ---
// (find_ring, dfs_find_cycle_generic, iupac_heterocyclo_name, best_heterocycle_numbering,
// find_longest_carbon_chain, dfs_longest_path, identify_substituents, name_substituent,
// count_branch_carbons, choose_numbering, assign_substituent_positions, format_substituents,
// is_aldehyde, find_carbon_ring, dfs_find_cycle, iupac_cyclo_name, identify_substituents_on_ring,
// best_ring_numbering, assign_substituents_on_ring, detect_ester, detect_ester_global)
// ... (keep these as in your current code) ...

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::UnGraph;

    // Test 1: Methane – a single carbon atom.
    #[test]
    fn test_methane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        graph.add_node(Element::C);
        let name = iupac_name(&graph);
        assert_eq!(name, "methane");
    }

    // Test 2: Ethane – two carbons with no substituents.
    #[test]
    fn test_ethane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        let a = graph.add_node(Element::C);
        let b = graph.add_node(Element::C);
        graph.add_edge(a, b, Bond::Single);
        let name = iupac_name(&graph);
        assert_eq!(name, "ethane");
    }

    // Test 3: 2-Methylpropane (isobutane) – a three-carbon chain with one methyl substituent on the middle carbon.
    #[test]
    fn test_2_methylpropane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Main chain: A - B - C
        let a = graph.add_node(Element::C);
        let b = graph.add_node(Element::C);
        let c = graph.add_node(Element::C);
        // Substituent (methyl) attached to B.
        let d = graph.add_node(Element::C);
        graph.add_edge(a, b, Bond::Single);
        graph.add_edge(b, c, Bond::Single);
        graph.add_edge(b, d, Bond::Single);

        let name = iupac_name(&graph);
        // Expected: "2-methyl-propane"
        assert_eq!(name, "2-methylpropane");
    }

    // Test 4: 2,3-Dimethylbutane – a four-carbon main chain with a methyl substituent on carbon 2 and carbon 3.
    #[test]
    fn test_2_3_dimethylbutane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Main chain: A - B - C - D
        let a = graph.add_node(Element::C);
        let b = graph.add_node(Element::C);
        let c = graph.add_node(Element::C);
        let d = graph.add_node(Element::C);
        // Substituents: methyl groups attached to B and C.
        let e = graph.add_node(Element::C);
        let f = graph.add_node(Element::C);
        graph.add_edge(a, b, Bond::Single);
        graph.add_edge(b, c, Bond::Single);
        graph.add_edge(c, d, Bond::Single);
        graph.add_edge(b, e, Bond::Single);
        graph.add_edge(c, f, Bond::Single);

        let name = iupac_name(&graph);
        // Expected: "2,3-dimethyl-butane"
        assert_eq!(name, "2,3-dimethylbutane");
    }

    // Test 5: 1-Chloroethane – a two-carbon chain with a chlorine substituent on carbon 1.
    #[test]
    fn test_1_chloroethane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        let a = graph.add_node(Element::C);
        let b = graph.add_node(Element::C);
        // Chlorine substituent attached to A.
        let c = graph.add_node(Element::Cl);
        graph.add_edge(a, b, Bond::Single);
        graph.add_edge(a, c, Bond::Single);

        let name = iupac_name(&graph);
        // Expected: "1-chloro-ethane"
        assert_eq!(name, "1-chloroethane");
    }

    // Test 6: Structure that can be interpreted as either "2-ethyl-pentane" or "3-methyl-hexane".
    // According to IUPAC rules, the correct name is "3-methyl-hexane".
    #[test]
    fn test_3_methyl_hexane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Construct a molecule that can be seen as having either a five-carbon chain with an ethyl group
        // or a six-carbon chain with a methyl group.
        // Here we construct it so that the algorithm selects a six-carbon main chain.
        //
        // Main chain: A - B - C - D - E - F
        let a = graph.add_node(Element::C);
        let b = graph.add_node(Element::C);
        let c = graph.add_node(Element::C);
        let d = graph.add_node(Element::C);
        let e = graph.add_node(Element::C);
        let f = graph.add_node(Element::C);
        // Branch: a methyl substituent attached to D.
        let g = graph.add_node(Element::C);
        graph.add_edge(a, b, Bond::Single);
        graph.add_edge(b, c, Bond::Single);
        graph.add_edge(c, d, Bond::Single);
        graph.add_edge(d, e, Bond::Single);
        graph.add_edge(e, f, Bond::Single);
        graph.add_edge(d, g, Bond::Single);

        let name = iupac_name(&graph);
        // Expected: "3-methyl-hexane"
        assert_eq!(name, "3-methylhexane");
    }

    // Test cyclohexane: 6 carbons in a ring with no substituents.
    #[test]
    fn test_cyclohexane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Create 6 carbon nodes.
        let nodes: Vec<_> = (0..6).map(|_| graph.add_node(Element::C)).collect();
        // Connect them in a cycle.
        for i in 0..6 {
            let next = (i + 1) % 6;
            graph.add_edge(nodes[i], nodes[next], Bond::Single);
        }
        let name = iupac_name(&graph);
        assert_eq!(name, "cyclohexane");
    }

    // Test methylcyclohexane: cyclohexane ring with one methyl substituent.
    #[test]
    fn test_methylcyclohexane() {
        init_logging("trace");
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Create cyclohexane ring nodes.
        let ring_nodes: Vec<_> = (0..6).map(|_| graph.add_node(Element::C)).collect();
        for i in 0..6 {
            let next = (i + 1) % 6;
            graph.add_edge(ring_nodes[i], ring_nodes[next], Bond::Single);
        }
        // Attach a methyl substituent (a single carbon) to one ring carbon.
        let methyl = graph.add_node(Element::C);
        graph.add_edge(ring_nodes[0], methyl, Bond::Single);
        let name = iupac_name(&graph);
        // Expected: "1-methyl-cyclohexane"
        assert_eq!(name, "1-methylcyclohexane");
    }

    // Test cyclopentanol: cyclopentane ring with a hydroxy substituent.
    #[test]
    fn test_cyclopentanol() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Create cyclopentane ring nodes (5 carbons).
        let ring_nodes: Vec<_> = (0..5).map(|_| graph.add_node(Element::C)).collect();
        for i in 0..5 {
            let next = (i + 1) % 5;
            graph.add_edge(ring_nodes[i], ring_nodes[next], Bond::Single);
        }
        // Attach an OH group (oxygen) to one ring carbon.
        let hydroxy = graph.add_node(Element::O);
        graph.add_edge(ring_nodes[0], hydroxy, Bond::Single);
        let name = iupac_name(&graph);
        // With symmetry, the substituent should get the lowest number: "1-hydroxy-cyclopentane".
        assert_eq!(name, "1-hydroxycyclopentane");
    }

    // Test 1,2-dimethylcyclopentane: cyclopentane with methyl substituents on adjacent carbons.
    #[test]
    fn test_dimethylcyclopentane() {
        let mut graph: MoleculeGraph = UnGraph::new_undirected();
        // Create cyclopentane ring nodes.
        let ring_nodes: Vec<_> = (0..5).map(|_| graph.add_node(Element::C)).collect();
        for i in 0..5 {
            let next = (i + 1) % 5;
            graph.add_edge(ring_nodes[i], ring_nodes[next], Bond::Single);
        }
        // Attach methyl substituents to two adjacent ring carbons.
        let methyl1 = graph.add_node(Element::C);
        let methyl2 = graph.add_node(Element::C);
        graph.add_edge(ring_nodes[0], methyl1, Bond::Single);
        graph.add_edge(ring_nodes[1], methyl2, Bond::Single);
        let name = iupac_name(&graph);
        // Expected: "1,2-dimethyl-cyclopentane"
        assert_eq!(name, "1,2-dimethylcyclopentane");
    }

    #[test]
    fn test_smiles_to_name() {
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            ("C", "methane"),
            ("CC", "ethane"),
            ("CCC", "propane"),
            ("CC(C)C", "2-methylpropane"),
            ("CC(C)C(C)C", "2,3-dimethylbutane"),
            ("CC(CC)CCC", "3-methylhexane"),
            ("C1CCCCC1", "cyclohexane"),
            ("C1CC(C)CCC1", "1-methylcyclohexane"),
            // Examples with oxygen:
            ("CCO", "1-hydroxyethane"), // Ethanol represented as hydroxy-ethane
            // Examples with nitrogen:
            ("CN", "1-aminomethane"), // Methylamine represented as amino-methane
            ("CCN", "1-aminoethane"), // Ethylamine represented as amino-ethane
            // Combined oxygen and nitrogen:
            ("CC(N)O", "1-amino-1-hydroxyethane"), // A 2-carbon chain with both NH2 and OH on the same carbon
            ("C=C", "ethene"),
            ("C#C", "ethyne"),
            ("C=CC=C", "1,3-butadiene"),
            ("c1ccccc1", "benzene"),
            ("Cc1ccccc1", "1-methylbenzene"),
            ("CC#N", "1-cyanoethane"),
            // Esters
            // Methyl acetate: CH₃COOCH₃
            ("COC(C)=O", "methyl-ethanoate"),
            // Ethyl acetate: CH₃COOC₂H₅
            ("CC(=O)OCC", "ethyl-ethanoate"),
            // Propyl propanoate: C₂H₅COOC₃H₇
            ("CCC(=O)OCCC", "propyl-propanoate"),
            // Ethers
            // Dimethyl ether: CH₃OCH₃
            ("COC", "1-methoxy-methane"),
            // Aldehydes
            // Formaldehyde: H–C(=O)–H
            ("C=O", "methanal"),
            // Acetaldehyde: CH₃CHO
            ("CC=O", "ethanal"),
            // Propionaldehyde: CH₃CH₂CHO
            ("C(CC)=O", "propanal"),
            // Isobutyraldehyde: (CH₃)₂CHCHO
            ("CC(C)C=O", "2-methylpropanal"),
            ("CC(O)C", "2-hydroxypropane"),
        ];

        for (smiles, correct_name) in smiles_and_correct_names {
            // Parse the SMILES string into a molecule graph.
            let molecule = parse_smiles(smiles)
                .unwrap_or_else(|_| panic!("Failed to parse SMILES: {}", smiles));

            // Optionally visualize the molecule graph.
            let mut dot_filename = String::new();
            write!(&mut dot_filename, "{}.dot", correct_name).unwrap();
            let mut png_filename = String::new();
            write!(&mut png_filename, "{}.png", correct_name).unwrap();
            visualize_graph(&molecule, &dot_filename, Some(&png_filename))
                .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: {}", smiles));

            // Generate the IUPAC name for the molecule.
            let generated_name = iupac_name(&molecule);
            assert_eq!(
                generated_name, correct_name,
                "For SMILES '{}' expected '{}' but got '{}'",
                smiles, correct_name, generated_name
            );
        }
    }

    #[test]
    fn test_esters() {
        // init_logging("trace");
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            ("COC(C)=O", "methyl-ethanoate"),
            ("CC(=O)OCC", "ethyl-ethanoate"),
            ("CCC(=O)OCCC", "propyl-propanoate"),
            ("CCCCOC(=O)C", "butyl-ethanoate"),
            ("C(CC)(=O)OC", "methyl-propanoate"),
            ("C(CC)(=O)OCC", "ethyl-propanoate"),
            ("C(CCC)(=O)OC", "methyl-butanoate"),
            ("C(CCC)(=O)OCC", "ethyl-butanoate"),
            ("O=C(OC)c1ccccc1", "methyl-benzoate"),
            ("O=C(OCC)c1ccccc1", "ethyl-benzoate"),
            ("CC(C)(C)OC(=O)C", "tert-butyl-ethanoate"),
            ("O=C(OC1=CC=CC=C1)C", "phenyl-ethanoate"),
            ("COC(=O)C1=CC=CO1", "methyl-furan-2-carboxylate"),
            ("COC(=O)C1=CC=CS1", "methyl-thiophene-2-carboxylate"),
            ("COC(=O)C1=NC=CC=C1", "methyl-pyridine-2-carboxylate"),
            ("C(C)OC(=O)C=1C=NC=CC1", "ethyl-pyridine-3-carboxylate"),
            ("O=C(OCCC)c1ccccc1", "propyl-benzoate"),
            ("C(C)(C)OC(C1=CC=CC=C1)=O", "isopropyl-benzoate"),
            ("COC(C(CCCC)CC)=O", "methyl-2-ethylhexanoate"),
            ("O=C(OC)C(=O)OC", "dimethyl-oxalate"),
        ];

        for (smiles, correct_name) in smiles_and_correct_names {
            // Parse the SMILES string into a molecule graph.
            let molecule = parse_smiles(smiles)
                .unwrap_or_else(|_| panic!("Failed to parse SMILES: {}", smiles));

            // Optionally visualize the molecule graph.
            let mut dot_filename = String::new();
            write!(&mut dot_filename, "{}.dot", correct_name).unwrap();
            let mut png_filename = String::new();
            write!(&mut png_filename, "{}.png", correct_name).unwrap();
            visualize_graph(&molecule, &dot_filename, Some(&png_filename))
                .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: {}", smiles));

            // Generate the IUPAC name for the molecule.
            let generated_name = iupac_name(&molecule);
            assert_eq!(
                generated_name, correct_name,
                "For SMILES '{}' expected '{}' but got '{}' in {:#?}",
                smiles, correct_name, generated_name, molecule
            );
            println!("Correctly named {:<26} as {:<32} ✅", smiles, generated_name);
        }
    }

    #[test]
    fn test_ethers() {
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            // Ethers
            ("COC", "1-methoxy-methane"),
            ("C(C)OCC", "1-ethoxy-ethane"),
            ("O1CCOCC1", "1,4-dioxacyclohexane"),
        ];

        for (smiles, correct_name) in smiles_and_correct_names {
            // Parse the SMILES string into a molecule graph.
            let molecule = parse_smiles(smiles)
                .unwrap_or_else(|_| panic!("Failed to parse SMILES: {}", smiles));

            // Optionally visualize the molecule graph.
            let mut dot_filename = String::new();
            write!(&mut dot_filename, "{}.dot", correct_name).unwrap();
            let mut png_filename = String::new();
            write!(&mut png_filename, "{}.png", correct_name).unwrap();
            visualize_graph(&molecule, &dot_filename, Some(&png_filename))
                .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: {}", smiles));

            // Generate the IUPAC name for the molecule.
            let generated_name = iupac_name(&molecule);
            assert_eq!(
                generated_name, correct_name,
                "For SMILES '{}' expected '{}' but got '{}'",
                smiles, correct_name, generated_name
            );
        }
    }

    #[test]
    fn test_aldehydes() {
        // init_logging(tracing::metadata::LevelFilter::TRACE);
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            // Aldehydes
            ("C=O", "methanal"),
            ("CC=O", "ethanal"),
            ("CCC=O", "propanal"),
            ("CC(C)C=O", "2-methylpropanal"),
            ("O1C(=CC=C1)C=O", "furan-2-carboxaldehyde"),
        ];

        for (smiles, correct_name) in smiles_and_correct_names {
            // Parse the SMILES string into a molecule graph.
            let molecule = parse_smiles(smiles)
                .unwrap_or_else(|_| panic!("Failed to parse SMILES: {}", smiles));

            // Optionally visualize the molecule graph.
            let mut dot_filename = String::new();
            write!(&mut dot_filename, "{}.dot", correct_name).unwrap();
            let mut png_filename = String::new();
            write!(&mut png_filename, "{}.png", correct_name).unwrap();
            visualize_graph(&molecule, &dot_filename, Some(&png_filename))
                .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: {}", smiles));

            // Generate the IUPAC name for the molecule.
            let generated_name = iupac_name(&molecule);
            assert_eq!(
                generated_name, correct_name,
                "For SMILES '{}' expected '{}' but got '{}'",
                smiles, correct_name, generated_name
            );
        }
    }

    #[test]
    fn test_pyrrolidine() {
        let molecule = parse_smiles("C1CCCN1")
            .unwrap_or_else(|_| panic!("Failed to parse SMILES: C1CCCN1"));
        let name = iupac_name(&molecule);
        assert_eq!(name, "pyrrolidine");
    }

    #[test]
    fn test_phenyl() {
        init_logging("trace");
        let molecule = parse_smiles("Cc1ccccc1")
            .unwrap_or_else(|_| panic!("Failed to parse SMILES: Cc1ccccc1"));
        visualize_graph(&molecule, "methylbenzene.dot", Some("methylbenzene.png"))
            .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: Cc1ccccc1"));
        let name = iupac_name(&molecule);
        assert_eq!(name, "1-methylbenzene");

        init_logging("trace");
        let molecule = parse_smiles("CCc1ccccc1")
            .unwrap_or_else(|_| panic!("Failed to parse SMILES: CCc1ccccc1"));
        visualize_graph(&molecule, "ethylbenzene.dot", Some("ethylbenzene.png"))
            .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: CCc1ccccc1"));
        let name = iupac_name(&molecule);
        assert_eq!(name, "1-ethylbenzene");
        
        init_logging("trace");
        let molecule = parse_smiles("C1(=CC=CC=C1)C1=CC=CC=C1")
            .unwrap_or_else(|_| panic!("Failed to parse SMILES: C1(=CC=CC=C1)C1=CC=CC=C1"));
        visualize_graph(&molecule, "phenylethane.dot", Some("phenylethane.png"))
            .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: C1(=CC=CC=C1)C1=CC=CC=C1"));
        let name = iupac_name(&molecule);
        assert_eq!(name, "1-phenylbenzene");
    }

    #[test]
    fn test_amino_acids() {
        // init_logging("trace");
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            ("NC(C(=O)O)CCCCN", "2,6-diaminohexanoic acid"),
            ("NC(C(=O)O)CCCNC(=N)N", "2-amino-5-guanidinopentanoic acid"),
            ("NC(C(=O)O)C", "2-aminopropanoic acid"),
            ("NC(C(=O)O)CC(N)=O", "2-amino-3-carbamoylpropanoic acid"),
            ("NC(C(=O)O)CC(=O)O", "2-aminobutanedioic acid"),
            ("NC(C(=O)O)CS", "2-amino-3-sulfhydrylpropanoic acid"),
            ("O=C(N)CCC(N)C(=O)O", "2,5-diamino-5-oxopentanoic acid"),
            ("NC(C(=O)O)CCC(=O)O", "2-aminopentanedioic acid"),
            ("C(C(=O)O)N", "2-aminoethanoic acid"),
            ("NC(C(=O)O)C(CC)C", "2-amino-3-methylpentanoic acid"),
            ("N1C(CCC1)C(=O)O", "pyrrolidine-2-carboxylic acid"),
            ("NC(C(=O)O)CO", "2-amino-3-hydroxypropanoic acid"),
            ("NC(C(=O)O)C(C)O", "2-amino-3-hydroxybutanoic acid"),
            ("NC(C(=O)O)C(C)C", "2-amino-3-methylbutanoic acid"),
        ];


        // ("NC(C(=O)O)CC1=CC=CC=C1", "2-amino-3-phenylpropanoic acid"),
        for (smiles, correct_name) in smiles_and_correct_names {
            // Parse the SMILES string into a molecule graph.
            let molecule = parse_smiles(smiles)
                .unwrap_or_else(|_| panic!("Failed to parse SMILES: {}", smiles));

            // Optionally visualize the molecule graph.
            let mut dot_filename = String::new();
            write!(&mut dot_filename, "{}.dot", correct_name).unwrap();
            let mut png_filename = String::new();
            write!(&mut png_filename, "{}.png", correct_name).unwrap();
            visualize_graph(&molecule, &dot_filename, Some(&png_filename))
                .unwrap_or_else(|_| panic!("Failed to visualize graph for SMILES: {}", smiles));

            // Generate the IUPAC name for the molecule.
            let generated_name = iupac_name(&molecule);
            assert_eq!(
                generated_name, correct_name,
                "For SMILES '{}' expected '{}' but got '{}'",
                smiles, correct_name, generated_name
            );
            println!("Correctly named {:<20} as {:<35} ✅", smiles, generated_name);
        }
    }
}
