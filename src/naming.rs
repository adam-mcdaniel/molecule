use petgraph::graph::{NodeIndex, UnGraph};
use std::collections::{BTreeMap, HashSet, HashMap};
use super::*;
use petgraph::visit::EdgeRef;
use tracing::debug;  // added logging from tracing

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
                            },
                            Some(Bond::Single) => {
                                // Check if oxygen connects to a carbon (other than node)
                                for nbr2 in graph.neighbors(nbr) {
                                    if nbr2 != node && (graph[nbr2].is_carbon()) {
                                        alkoxy_oxygen = Some(nbr);
                                        debug!("Found alkoxy oxygen {:?} connected to node {:?} from {:?}", nbr, node, nbr2);
                                        break;
                                    }
                                }
                            },
                            _ => {}
                        }
                    }
                }
            }
            if has_carbonyl && alkoxy_oxygen.is_some() {
                debug!("Ester detected: acyl carbon {:?} with alkoxy oxygen {:?}", node, alkoxy_oxygen.unwrap());
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
fn detect_ether_global(graph: &MoleculeGraph) -> Option<(Vec<NodeIndex>, Vec<NodeIndex>, NodeIndex)> {
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
                        debug!("Oxygen node {:?} is part of a carbonyl with neighbor {:?}", o, nbr);
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
            debug!("Ether candidate at oxygen {:?}: comp1 size {}, comp2 size {}", o, comp1.len(), comp2.len());
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

/// Names an ester molecule using the globally detected ester atoms.
fn iupac_ester_name(graph: &MoleculeGraph, acyl_carbon: NodeIndex, alkoxy_oxygen: NodeIndex) -> String {
    debug!("Naming ester: acyl carbon {:?}, alkoxy oxygen {:?}", acyl_carbon, alkoxy_oxygen);
    let acyl_name = name_acyl_group(graph, acyl_carbon);
    if let Some(alkoxy_carb) = graph.neighbors(alkoxy_oxygen)
        .find(|&nbr| nbr != acyl_carbon && (graph[nbr].is_carbon()))
    {
        if let Some(ring_name) = alkoxy_ring_name(graph, alkoxy_carb) {
            debug!("Using aromatic alkoxy ring name: {}", ring_name);
            return format!("{}-{}", ring_name, acyl_name);
        }
        if is_tert_butyl(graph, alkoxy_carb, alkoxy_oxygen) {
            debug!("Detected tert-butyl substituent");
            return format!("tert-butyl-{}", acyl_name);
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
        debug!("Ester naming: using alkoxy prefix {}", alkoxy_name);
        return format!("{}-{}", alkoxy_name, acyl_name);
    }
    acyl_name
}

/// Names an aldehyde molecule given its longest carbon chain and the index
/// (in that chain) at which the terminal aldehyde group occurs.
fn iupac_aldehyde_name(graph: &MoleculeGraph, mut main_chain: Vec<NodeIndex>, aldehyde_index: usize) -> String {
    debug!("Naming aldehyde: chain length {}, aldehyde index {}", main_chain.len(), aldehyde_index);
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
    debug!("Aldehyde naming: final chain length {}, substituent prefix '{}'", chain_length, substituent_prefix);
    format!("{}al", substituent_prefix + root)
}

/// Names an ether molecule given the two disconnected carbon fragments and the bridging oxygen.
/// The main_chain is chosen as the larger (or preferred) fragment.
fn iupac_ether_name(graph: &MoleculeGraph, mut main_chain: Vec<NodeIndex>, alkoxy_chain: Vec<NodeIndex>, oxygen: NodeIndex) -> String {
    debug!("Naming ether: main chain length {}, alkoxy chain length {}, oxygen {:?}", main_chain.len(), alkoxy_chain.len(), oxygen);
    if let Some(&attach) = main_chain.iter().find(|&&n| graph.find_edge(n, oxygen).is_some()) {
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

// --- Modified acyclic naming routine ---
pub fn iupac_acyclic_name(graph: &MoleculeGraph) -> String {
    debug!("Starting acyclic naming routine");
    if let Some((acyl_carbon, alkoxy_oxygen)) = detect_ester_global(graph) {
        debug!("Acyclic naming: Ester branch taken");
        return iupac_ester_name(graph, acyl_carbon, alkoxy_oxygen);
    }
    if let Some((main_chain, aldehyde_index)) = detect_aldehyde_global(graph) {
        debug!("Acyclic naming: Aldehyde branch taken");
        return iupac_aldehyde_name(graph, main_chain, aldehyde_index);
    }
    if let Some((main_chain, alkoxy_chain, oxygen)) = detect_ether_global(graph) {
        debug!("Acyclic naming: Ether branch taken");
        return iupac_ether_name(graph, main_chain, alkoxy_chain, oxygen);
    }
    // NEW: First check for an exocyclic aldehyde attached to a heterocycle.
    if let Some(name) = detect_exocyclic_aldehyde(graph) {
        debug!("Acyclic naming: Exocyclic aldehyde branch taken");
        return name;
    }
    
    debug!("Acyclic naming: Falling back to standard acyclic naming");
    let main_chain = find_longest_carbon_chain(graph);
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
        },
    };
    let subs_final = assign_substituent_positions(&final_main_chain, &subs_original, numbering_order);
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
    
    // Determine unsaturation (double and triple bonds) along the main chain.
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
                },
                Some(Bond::Triple) => {
                    if chain_length > 2 {
                        unsat_locs.push(i + 1);
                    }
                    triple_count += 1;
                },
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
    debug!("Acyclic naming: final name composed with prefix '{}' and suffix '{}'", substituent_prefix, suffix);
    format!("{}{}{}{}", substituent_prefix, unsat_prefix, unsaturated_base, suffix)
}

/// --- Main naming routine ---
pub fn iupac_name(graph: &MoleculeGraph) -> String {
    debug!("Starting IUPAC naming process");
    if detect_ester_global(graph).is_some() {
        debug!("Ester detected globally; using acyclic naming branch");
        return iupac_acyclic_name(graph);
    }
    // NEW: Prioritize exocyclic aldehyde detection if present.
    if let Some(name) = detect_exocyclic_aldehyde(graph) {
        debug!("Exocyclic aldehyde detected; using that naming branch");
        return name;
    }
    if let Some(ring) = find_ring(graph) {
        debug!("Ring detected with {} atoms", ring.len());
        if ring.len() >= 3 {
            if ring.iter().all(|&n| graph[n].is_carbon()) {
                debug!("Naming as cyclic alkane");
                return iupac_cyclo_name(graph, &ring);
            } else {
                debug!("Naming as heterocycle");
                return iupac_heterocyclo_name(graph, &ring);
            }
        }
    }
    debug!("No cyclic structure detected; falling back to acyclic naming");
    iupac_acyclic_name(graph)
}

// --- NEW: Generic ring detection that does not filter out heteroatoms ---
fn find_ring(graph: &MoleculeGraph) -> Option<Vec<NodeIndex>> {
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    for node in graph.node_indices() {
        if !visited.contains(&node) {
            if let Some(cycle) = dfs_find_cycle_generic(graph, node, None, &mut stack, &mut visited) {
                return Some(cycle);
            }
        }
    }
    None
}

fn dfs_find_cycle_generic(
    graph: &MoleculeGraph,
    current: NodeIndex,
    parent: Option<NodeIndex>,
    stack: &mut Vec<NodeIndex>,
    visited: &mut HashSet<NodeIndex>,
) -> Option<Vec<NodeIndex>> {
    visited.insert(current);
    stack.push(current);
    for neighbor in graph.neighbors(current) {
        if Some(neighbor) == parent {
            continue;
        }
        if stack.contains(&neighbor) {
            let pos = stack.iter().position(|&x| x == neighbor).unwrap();
            return Some(stack[pos..].to_vec());
        } else if !visited.contains(&neighbor) {
            if let Some(cycle) = dfs_find_cycle_generic(graph, neighbor, Some(current), stack, visited) {
                return Some(cycle);
            }
        }
    }
    stack.pop();
    None
}

// --- NEW: Heterocyclic naming function ---
fn iupac_heterocyclo_name(graph: &MoleculeGraph, ring: &Vec<NodeIndex>) -> String {
    let n = ring.len();
    // Determine the parent cycloalkane name based on ring size.
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

    // Choose the best numbering for the ring so that heteroatoms get low locants.
    let (_best_order, hetero_atoms) = best_heterocycle_numbering(ring, graph);
    // Group heteroatoms by type.
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
    // For simplicity, if there is only one type (as in our test case)
    // build a prefix such as "1,4-dioxa" (for two oxygens).
    let prefix = if groups.len() == 1 {
        let (group, mut locants) = groups.into_iter().next().unwrap();
        locants.sort();
        let locant_str = locants.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
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
    let mut best_locants: Vec<usize> = ring.iter().enumerate()
        .filter(|&(_, &node)| !graph[node].is_carbon())
        .map(|(i, _)| i + 1)
        .collect();
    best_locants.sort();

    let mut best_hetero: Vec<(usize, Element)> = ring.iter().enumerate()
         .filter(|&(_, &node)| !graph[node].is_carbon())
         .map(|(i, &node)| (i + 1, graph[node])).collect();
    best_hetero.sort_by_key(|&(loc, _)| loc);

    // Try all rotations and directions.
    for start in 0..n {
        let candidate: Vec<NodeIndex> = ring[start..].iter()
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
fn identify_substituents(
    graph: &MoleculeGraph,
    main_chain: &[NodeIndex],
) -> Vec<(usize, String)> {
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

/// Adjusted naming for substituents:
/// - A double-bonded oxygen (to the parent) is named "oxo" (for aldehydes).
/// - A single-bonded oxygen that connects to another carbon is named as an ether (e.g. "methoxy").
/// - Otherwise, oxygen is "hydroxy".
fn name_substituent(graph: &MoleculeGraph, neighbor: NodeIndex, exclude: NodeIndex) -> String {
    match graph[neighbor] {
        Element::C => {
            // For carbon branches, count how many connected carbon atoms are in this substituent.
            let branch_size = count_branch_carbons(graph, neighbor, exclude);
            match branch_size {
                1 => "methyl".to_string(),
                2 => "ethyl".to_string(),
                3 => "propyl".to_string(),
                4 => "butyl".to_string(),
                5 => "pentyl".to_string(),
                _ => format!("{}-carbon alkyl", branch_size),
            }
        },
        Element::O => {
            // Check the bond type between the parent and oxygen.
            if let Some(edge) = graph.find_edge(exclude, neighbor) {
                if graph.edge_weight(edge) == Some(&Bond::Double) {
                    // A double bond indicates a carbonyl (aldehyde) group.
                    return "oxo".to_string();
                }
            }
            // For single-bonded oxygen, determine if it connects to another carbon.
            let carbon_neighbors: Vec<_> = graph.neighbors(neighbor)
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
        },
        Element::N => {
            if let Some(edge) = graph.find_edge(exclude, neighbor) {
                if let Some(bond) = graph.edge_weight(edge) {
                    if *bond == Bond::Triple {
                        return "cyano".to_string();
                    }
                }
            }
            "amino".to_string()
        },
        Element::F => "fluoro".to_string(),
        Element::Cl => "chloro".to_string(),
        Element::Br => "bromo".to_string(),
        _ => "substituent".to_string(),
    }
}

fn count_branch_carbons(
    graph: &MoleculeGraph,
    start: NodeIndex,
    exclude: NodeIndex,
) -> usize {
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

fn choose_numbering(
    main_chain: &[NodeIndex],
    substituents: &[(usize, String)],
) -> NumberingOrder {
    if substituents.is_empty() {
        return NumberingOrder::Original;
    }
    let n = main_chain.len();
    let mut forward_locants: Vec<usize> = substituents.iter().map(|(pos, _)| *pos).collect();
    forward_locants.sort();
    let mut reverse_locants: Vec<usize> = substituents.iter().map(|(pos, _)| n - pos + 1).collect();
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
                NumberingOrder::Reversed => n - pos + 1,
            };
            Substituent {
                position: new_pos,
                name: name.clone(),
            }
        })
        .collect()
}

// Updated substituent formatting function.
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
        parts.push(format!("{}-{}{}", locant_str, prefix, name));
    }
    // Join groups with commas and append a trailing hyphen.
    format!("{}-", parts.join(","))
}


fn find_carbon_ring(graph: &MoleculeGraph) -> Option<Vec<NodeIndex>> {
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    for node in graph.node_indices() {
        if graph[node].is_carbon() && !visited.contains(&node) {
            if let Some(cycle) = dfs_find_cycle(graph, node, None, &mut stack, &mut visited) {
                return Some(cycle);
            }
        }
    }
    None
}

fn dfs_find_cycle(
    graph: &MoleculeGraph,
    current: NodeIndex,
    parent: Option<NodeIndex>,
    stack: &mut Vec<NodeIndex>,
    visited: &mut HashSet<NodeIndex>,
) -> Option<Vec<NodeIndex>> {
    visited.insert(current);
    stack.push(current);
    for neighbor in graph.neighbors(current) {
        // Only consider carbon atoms (including aromatic carbons) for the ring.
        if !graph[neighbor].is_carbon() {
            continue;
        }
        if Some(neighbor) == parent {
            continue;
        }
        if stack.contains(&neighbor) {
            // Found a cycle: extract the portion of the stack from neighbor to current.
            let pos = stack.iter().position(|&x| x == neighbor).unwrap();
            return Some(stack[pos..].to_vec());
        } else if !visited.contains(&neighbor) {
            if let Some(cycle) = dfs_find_cycle(graph, neighbor, Some(current), stack, visited) {
                return Some(cycle);
            }
        }
    }
    stack.pop();
    None
}

/// Returns the IUPAC name for a cyclic molecule (cycloalkane possibly with substituents).
fn iupac_cyclo_name(graph: &MoleculeGraph, ring: &Vec<NodeIndex>) -> String {
    // Check if the ring is fully aromatic.
    let is_aromatic_ring = ring.iter().all(|&n| graph[n].is_aromatic() && graph[n].is_carbon());
    
    if is_aromatic_ring && ring.len() == 6 {
        // Identify substituents on the aromatic ring.
        let subs = identify_substituents_on_ring(graph, ring);
        if subs.is_empty() {
            return "benzene".to_string();
        } else {
            // Determine the best numbering for the ring by trying all rotations and directions.
            let (best_ring, best_subs) = best_ring_numbering(ring, graph, &subs);
            let substituent_prefix = format_substituents(&best_subs);
            return format!("{}benzene", substituent_prefix);
        }
    }
    
    // Non-aromatic (or aromatic but not 6-membered) cyclic naming.
    let subs = identify_substituents_on_ring(graph, ring);
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
        let mut candidate = ring[start..].iter().chain(ring[..start].iter()).cloned().collect::<Vec<_>>();
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
            let mut candidate_locants: Vec<usize> = candidate_subs.iter().map(|s| s.position).collect();
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

// --- NEW: Acyl Chain Determination ---
// Starting at the acyl carbon, find the longest path that follows only carbon atoms.
// This chain represents the acid-derived portion of the ester.
fn find_acyl_chain(graph: &MoleculeGraph, acyl_carbon: NodeIndex) -> Vec<NodeIndex> {
    let mut visited = HashSet::new();
    dfs_longest_path(graph, acyl_carbon, &mut visited)
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
                    let subst = base[..base.len()-3].to_string() + "yl";
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
    let a = best_order.iter().position(|&node| node == attach).unwrap_or(0);
    let d1 = if a >= h { a - h } else { a + n - h };
    let d2 = n - d1;
    d1.min(d2) + 1
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
                if base != "benzene" { // only treat heterocycles specially
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
                if base != "benzene" { // only for heterocycles
                    let (best_order, _) = best_heterocycle_numbering(&ring, graph);
                    let locant = conventional_locant_for_heterocycle(graph, &best_order, nbr);
                    return Some(format!("{}-{}-carboxaldehyde", base, locant));
                }
            }
        }
    }
    None
}

/// Modified general acyl group naming function for esters.
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
        return name_heterocyclic_acid(graph, &ring, acyl_carbon);
    }

    // Otherwise, check if one of its non-oxygen neighbors is in a heterocyclic ring.
    if let Some(exo_name) = name_exocyclic_acid(graph, acyl_carbon) {
        return exo_name;
    }

    if let Some(name) = name_exocyclic_aldehyde(graph, acyl_carbon) {
        return name;
    }

    // Fallback: treat the acyl group as acyclic.
    let acyl_chain = find_acyl_chain(graph, acyl_carbon);
    let acyl_length = acyl_chain.len();
    let acyl_base = match acyl_length {
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
        n => return format!("{}-carbon acid chain (not supported)", n),
    };
    format!("{}oate", acyl_base)
}

/// --- GENERAL HELPER FUNCTIONS FOR RING-BASED ACID NAMING ---

/// Returns a base name for a ring based on its elemental composition.
/// This is a data-driven lookup that minimizes hard coding.
/// Extend this map as needed.
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
    // For benzene: 6 carbons, no heteroatoms.
    if count_c == 6 && count_o == 0 && count_n == 0 && count_s == 0 {
        return Some("benzene".to_string());
    }
    // For pyridine: 5 carbons and one nitrogen.
    if count_c == 5 && count_n == 1 && count_o == 0 && count_s == 0 {
        return Some("pyridine".to_string());
    }
    // For furan: 4 carbons and one oxygen.
    if count_c == 4 && count_o == 1 && count_n == 0 && count_s == 0 {
        return Some("furan".to_string());
    }
    // For thiophene: 4 carbons and one sulfur.
    if count_c == 4 && count_s == 1 && count_o == 0 && count_n == 0 {
        return Some("thiophene".to_string());
    }
    // Fallback: build a generic name by concatenating element symbols and counts.
    let mut freq: HashMap<String, usize> = HashMap::new();
    for &n in ring {
        let sym = graph[n].symbol().to_string();
        *freq.entry(sym).or_insert(0) += 1;
    }
    let mut parts: Vec<String> = freq.into_iter().map(|(s, c)| format!("{}{}", s, c)).collect();
    parts.sort();
    Some(parts.join(""))
}

/// Given a ring that contains the acyl carbon, returns a systematic acid name
/// of the form "<base>-<locant>-carboxylate". The numbering is determined by the
/// existing best_heterocycle_numbering routine.
fn name_heterocyclic_acid(graph: &MoleculeGraph, ring: &[NodeIndex], acyl_carbon: NodeIndex) -> String {
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
    if let Some(ring) = find_ring(graph) {
        if ring.contains(&node) {
            return Some(ring);
        }
    }
    None
}

/// --- Alkoxy Naming Helpers ---
/// Detects whether the candidate alkoxy carbon is a tert-butyl group.
fn is_tert_butyl(graph: &MoleculeGraph, node: NodeIndex, exclude: NodeIndex) -> bool {
    if !graph[node].is_carbon() {
        return false;
    }
    let carbon_neighbors: Vec<NodeIndex> = graph.neighbors(node)
        .filter(|&nbr| nbr != exclude && graph[nbr].is_carbon())
        .collect();
    if carbon_neighbors.len() != 3 {
        return false;
    }
    // Check that each neighbor is terminal.
    for nbr in carbon_neighbors {
        let count = graph.neighbors(nbr)
            .filter(|&x| x != node && graph[x].is_carbon())
            .count();
        if count != 0 {
            return false;
        }
    }
    true
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
        assert_eq!(name, "2-methyl-propane");
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
        assert_eq!(name, "2,3-dimethyl-butane");
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
        assert_eq!(name, "1-chloro-ethane");
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
        assert_eq!(name, "3-methyl-hexane");
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
        assert_eq!(name, "1-methyl-cyclohexane");
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
        assert_eq!(name, "1-hydroxy-cyclopentane");
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
        assert_eq!(name, "1,2-dimethyl-cyclopentane");
    }

    #[test]
    fn test_smiles_to_name() {
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            ("C", "methane"),
            ("CC", "ethane"),
            ("CCC", "propane"),
            ("CC(C)C", "2-methyl-propane"),
            ("CC(C)C(C)C", "2,3-dimethyl-butane"),
            ("CC(CC)CCC", "3-methyl-hexane"),
            ("C1CCCCC1", "cyclohexane"),
            ("C1CC(C)CCC1", "1-methyl-cyclohexane"),
            // Examples with oxygen:
            ("CCO", "1-hydroxy-ethane"),         // Ethanol represented as hydroxy-ethane
            // Examples with nitrogen:
            ("CN", "1-amino-methane"),           // Methylamine represented as amino-methane
            ("CCN", "1-amino-ethane"),           // Ethylamine represented as amino-ethane
            // Combined oxygen and nitrogen:
            ("CC(N)O", "1-amino,1-hydroxy-ethane"), // A 2-carbon chain with both NH2 and OH on the same carbon

            ("C=C", "ethene"),
            ("C#C", "ethyne"),
            ("C=CC=C", "1,3-butadiene"),
            ("c1ccccc1", "benzene"),
            ("Cc1ccccc1", "1-methyl-benzene"),
            ("CC#N", "1-cyano-ethane"),

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
            ("CC(C)C=O", "2-methyl-propanal"),
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
            ("O=C(OC)c1ccccc1","methyl-benzoate"),
            ("O=C(OCC)c1ccccc1", "ethyl-benzoate"),
            ("CC(C)(C)OC(=O)C", "tert-butyl-ethanoate"),
            ("O=C(OC1=CC=CC=C1)C", "phenyl-ethanoate"),
            ("COC(=O)C1=CC=CO1", "methyl-furan-2-carboxylate"),
            ("COC(=O)C1=CC=CS1", "methyl-thiophene-2-carboxylate"),
            ("COC(=O)C1=NC=CC=C1", "methyl-pyridine-2-carboxylate"),
            ("C(C)OC(=O)C=1C=NC=CC1", "ethyl-pyridine-3-carboxylate"),
            ("O=C(OCCC)c1ccccc1", "propyl-benzoate"),
            // ("C(C)(C)OC(C1=CC=CC=C1)=O", "isopropyl-benzoate"),
            // ("COC(C(CCCC)CC)=O", "methyl-2-ethylhexanoate")
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
    fn test_ethers() {
        // A list of tuples: (SMILES string, expected IUPAC name)
        let smiles_and_correct_names = vec![
            // Ethers
            ("COC", "1-methoxy-methane"),
            ("C(C)OCC", "1-ethoxy-ethane"),
            ("O1CCOCC1", "1,4-dioxacyclohexane")
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
            ("CC(C)C=O", "2-methyl-propanal"),
            ("O1C(=CC=C1)C=O", "furan-2-carboxaldehyde")
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
}


