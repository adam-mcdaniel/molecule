mod smiles;
pub use smiles::*;

mod iupac;
pub use iupac::*;

/*
/// Returns a canonical labeling for the graph as a vector of NodeIndex in canonical order,
/// using a simplified partition-refinement/backtracking algorithm (in the spirit of Nauty).
/// This implementation is only practical for very small graphs.
pub fn nauty_canonical_labeling<T: Clone, E: Clone>(
    graph: &Graph<T, E, Undirected>
) -> Vec<NodeIndex> {
    let n = graph.node_count();

    // --- Initial Partition ---
    // We partition nodes by their degree.
    let mut groups: HashMap<usize, Vec<NodeIndex>> = HashMap::new();
    for node in graph.node_indices() {
        let deg = graph.neighbors(node).count();
        groups.entry(deg).or_default().push(node);
    }
    let mut partition: Vec<Vec<NodeIndex>> = groups.into_values().collect();
    // Sort each cell and then sort cells by their first element index.
    for cell in partition.iter_mut() {
        cell.sort_by_key(|&node| node.index());
    }
    partition.sort_by_key(|cell| cell[0].index());

    // --- Helper Functions ---

    /// Flatten a partition (Vec<Vec<NodeIndex>>) into a single Vec<NodeIndex>
    /// in the order given by the partition cells.
    fn flatten_partition(part: &Vec<Vec<NodeIndex>>) -> Vec<NodeIndex> {
        part.iter().flatten().copied().collect()
    }

    /// Given a permutation of the nodes (as a vector of NodeIndex),
    /// return a string representation of the graph's connectivity.
    /// Here we build an "upper-triangular" string of 0s and 1s.
    fn permutation_representation<T: Clone, E: Clone>(
        graph: &Graph<T, E, Undirected>,
        perm: &Vec<NodeIndex>
    ) -> String {
        let n = perm.len();
        let mut repr = String::new();
        for i in 0..n {
            for j in (i+1)..n {
                if graph.find_edge(perm[i], perm[j]).is_some() {
                    repr.push('1');
                } else {
                    repr.push('0');
                }
            }
        }
        repr
    }

    // We'll store the best (lexicographically smallest) labeling found.
    let mut best_labeling: Option<Vec<NodeIndex>> = None;
    let mut best_repr: Option<String> = None;

    /// Recursive search that refines the partition and branches on non-discrete cells.
    fn search<T: Clone, E: Clone>(
        graph: &Graph<T, E, Undirected>,
        partition: Vec<Vec<NodeIndex>>,
        best_labeling: &mut Option<Vec<NodeIndex>>,
        best_repr: &mut Option<String>
    ) {
        // If every cell is a singleton, we have a complete candidate labeling.
        if partition.iter().all(|cell| cell.len() == 1) {
            let candidate = flatten_partition(&partition);
            let repr = permutation_representation(graph, &candidate);
            if best_repr.is_none() || &repr < best_repr.as_ref().unwrap() {
                *best_repr = Some(repr);
                *best_labeling = Some(candidate);
            }
            return;
        }
        // Otherwise, choose the first non-discrete cell.
        let mut new_partition = partition.clone();
        let cell_idx = new_partition.iter().position(|cell| cell.len() > 1).unwrap();
        let cell = new_partition.remove(cell_idx);

        // Branch: for each element in the ambiguous cell, fix that element first.
        for &elem in &cell {
            let mut branch_partition = new_partition.clone();
            // Create one singleton cell for the chosen element.
            branch_partition.push(vec![elem]);
            // Create another cell for the remaining elements (if any).
            let remaining: Vec<NodeIndex> = cell.iter().filter(|&&x| x != elem).copied().collect();
            if !remaining.is_empty() {
                branch_partition.push(remaining);
            }
            // Sort cells to enforce an order.
            for cell in branch_partition.iter_mut() {
                cell.sort_by_key(|&n| n.index());
            }
            branch_partition.sort_by_key(|cell| cell[0].index());
            // Recursively search this branch.
            search(graph, branch_partition, best_labeling, best_repr);
        }
    }

    search(graph, partition, &mut best_labeling, &mut best_repr);
    best_labeling.unwrap_or_else(|| (0..n).map(|i| NodeIndex::new(i)).collect())
}

/// Converts a MoleculeGraph (with Element nodes and Bond edges) into a new Graph<(), (), Undirected>
/// that has the same topology (nodes and edges) but no attached data.
pub fn molecule_graph_to_unit_graph(mol: &MoleculeGraph) -> Graph<(), (), Undirected> {
    let mut new_graph = Graph::<(), (), Undirected>::new_undirected();
    // Map each node in the original graph to a node in the new graph.
    let mut node_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    
    // Iterate over all nodes in the original graph.
    for node in mol.node_indices() {
        let new_node = new_graph.add_node(());
        node_map.insert(node, new_node);
    }
    
    // Iterate over all edges in the original graph and add them to the new graph.
    for edge in mol.edge_references() {
        let source = node_map[&edge.source()];
        let target = node_map[&edge.target()];
        new_graph.add_edge(source, target, ());
    }
    
    new_graph
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_graph_to_unit_graph() {
        // let mol = test_molecule_graph();
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        let unit_graph = molecule_graph_to_unit_graph(&mol);
        
        // Check that the number of nodes is the same.
        assert_eq!(mol.node_count(), unit_graph.node_count());
        
        // Check that the number of edges is the same.
        assert_eq!(mol.edge_count(), unit_graph.edge_count());
        
        // Check that the topology is the same.
        for edge in mol.edge_references() {
            let source = edge.source();
            let target = edge.target();
            assert!(unit_graph.find_edge(source, target).is_some() || unit_graph.find_edge(target, source).is_some());
        }
    }

    #[test]
    fn test_canonize() {
        use nauty_pet::prelude::*;
        use nauty_pet::canon::*;
        let mol1 = parse_smiles("c1c(R)c1").unwrap();
        let mol2 = parse_smiles("c1cc1R").unwrap();
        
        // let canonized1 = mol1.into_canon();
        // let canonized2 = mol2.into_canon();
        // let canonized1 = mol1.clone();
        // let canonized2 = mol2.clone();
        let canonized1 = canonize(&mol1);
        let canonized2 = canonize(&mol2);
        
        println!("Canonized 1: {:?}", canonized1);
        println!("Canonized 2: {:?}", canonized2);

        // println!("Is identical: {}", canonized1.is_identical(&canonized2));

        let org1 = OrganicMolecule::from(canonized1);
        org1.visualize("canonized1.png").unwrap();

        let org2 = OrganicMolecule::from(canonized2);
        org2.visualize("canonized2.png").unwrap();


        // Convert to smiles
        let smiles1 = org1.to_smiles().unwrap();
        let smiles2 = org2.to_smiles().unwrap();

        println!("Smiles 1: {}", smiles1);
        println!("Smiles 2: {}", smiles2);
    }
}
*/
