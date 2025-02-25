use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use std::collections::{BTreeMap, HashSet};
use super::*;

/// Main entry point for naming a molecule. If a cyclic structure is detected,
/// the cyclic naming routine is used; otherwise, an acyclic (chain) name is produced.
pub fn iupac_name(graph: &MoleculeGraph) -> String {
    // Try to detect a ring consisting solely of carbons.
    if let Some(ring) = find_carbon_ring(graph) {
        // Use cyclic naming if the ring has at least 3 members.
        if ring.len() >= 3 {
            return iupac_cyclo_name(graph, &ring);
        }
    }
    // Fall back to acyclic naming.
    iupac_acyclic_name(graph)
}

/// === ACYCLIC NAMING (previous implementation) ===

pub fn iupac_acyclic_name(graph: &MoleculeGraph) -> String {
    let main_chain = find_longest_carbon_chain(graph);
    if main_chain.is_empty() {
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
    // Determine unsaturation along the main chain.
    let mut unsat_locs = Vec::new();
    let mut double_count = 0;
    let mut triple_count = 0;
    for i in 0..(final_main_chain.len() - 1) {
        let n1 = final_main_chain[i];
        let n2 = final_main_chain[i + 1];
        if let Some(edge) = graph.find_edge(n1, n2) {
            match graph.edge_weight(edge) {
                Some(Bond::Double) => {
                    // For chains longer than 2, record locants.
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

    // Choose the suffix based on unsaturation.
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

    // For chains with unsaturation and more than one double or triple bond, insert an "a" between the base and suffix.
    // (For example, for a 4-carbon chain with two double bonds: "but" becomes "buta", so that "buta" + "diene" becomes "butadiene".)
    let unsaturated_base = if (double_count > 1 || triple_count > 1) && chain_length >= 4 {
        format!("{}a", base_name)
    } else {
        base_name.to_string()
    };

    // If chain is two carbons or there are no unsaturations, omit the locant prefix.
    let unsat_prefix = if chain_length == 2 || (double_count == 0 && triple_count == 0) {
        "".to_string()
    } else {
        let locants: Vec<String> = unsat_locs.iter().map(|loc| loc.to_string()).collect();
        format!("{}-", locants.join(","))
    };

    format!("{}{}{}{}", substituent_prefix, unsat_prefix, unsaturated_base, suffix)
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

/// For each carbon in the main chain, identify any neighbor (not in the chain)
/// as a substituent.
fn identify_substituents(
    graph: &MoleculeGraph,
    main_chain: &Vec<NodeIndex>,
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

fn name_substituent(graph: &MoleculeGraph, neighbor: NodeIndex, exclude: NodeIndex) -> String {
    if graph[neighbor] == Element::C {
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
    } else {
        match graph[neighbor] {
            Element::O => "hydroxy".to_string(),
            Element::N => {
                // Check the bond type between the parent atom and this nitrogen.
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
                if !visited.contains(&neighbor) && graph[neighbor] == Element::C {
                    stack.push(neighbor);
                }
            }
        }
    }
    visited.len() - 1
}

fn choose_numbering(
    main_chain: &Vec<NodeIndex>,
    substituents: &Vec<(usize, String)>,
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
    let mut groups: std::collections::BTreeMap<String, Vec<usize>> = std::collections::BTreeMap::new();
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
    // Join groups with commas (instead of hyphens) and append a trailing hyphen.
    format!("{}-", parts.join(","))
}

/// === CYCLIC NAMING ===
/// 
// Helper function to check if an element is a carbon (either aliphatic or aromatic).
fn is_carbon(element: Element) -> bool {
    matches!(element, Element::C | Element::CAromatic)
}

fn find_carbon_ring(graph: &MoleculeGraph) -> Option<Vec<NodeIndex>> {
    let mut visited = HashSet::new();
    let mut stack = Vec::new();
    for node in graph.node_indices() {
        if is_carbon(graph[node]) && !visited.contains(&node) {
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
        if !is_carbon(graph[neighbor]) {
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
    let is_aromatic_ring = ring.iter().all(|&n| graph[n] == Element::CAromatic);
    
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

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::NodeIndex;
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