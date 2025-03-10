use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use petgraph::graph::{NodeIndex, UnGraph, };
use petgraph::visit::IntoNodeReferences;
use petgraph::visit::{EdgeRef, };
use petgraph::EdgeType;
use crate::*;

pub trait MorganCanonize {
    fn morgan_canonize(&self) -> Self;
}

impl MorganCanonize for MoleculeGraph {
    fn morgan_canonize(&self) -> Self {
        let labels = morgan_algorithm(self, 100);
        let canonical_graph = rebuild_canonical_graph(self, &labels);
        canonical_graph
    }
}

pub fn is_identical<N, E, T>(g1: &Graph<N, E, T>, g2: &Graph<N, E, T>) -> bool
where
    N: Ord + std::fmt::Debug,
    E: Ord + std::fmt::Debug,
    T: EdgeType,
{
    // Get all the node and edge indices
    let g1_nodes: Vec<_> = g1.node_indices().map(|n| &g2[n]).collect();
    let g2_nodes: Vec<_> = g2.node_indices().map(|n| &g2[n]).collect();

    let g1_edges: Vec<_> = g1.edge_references().map(|e| (e.source(), e.target())).collect();
    let g2_edges: Vec<_> = g2.edge_references().map(|e| (e.source(), e.target())).collect();

    // Check if the number of nodes and edges are the same
    if g1_nodes.len() != g2_nodes.len() || g1_edges.len() != g2_edges.len() {
        return false;
    }

    // Confirm all nodes are identical
    if g1_nodes != g2_nodes {
        return false;
    }

    if g1_edges != g2_edges {
        return false;
    }

    // Check if the node data is the same
    for (n1, n2) in g1_nodes.iter().zip(g2_nodes.iter()) {
        if n1 != n2 {
            return false;
        }
    }

    // Check if the edge data is the same
    for (e1, e2) in g1_edges.iter().zip(g2_edges.iter()) {
        if e1 != e2 {
            return false;
        }
    }

    true
}
/// Computes a hash value for any hashable object.
fn compute_hash<T: Hash>(t: &T) -> usize {
    let mut hasher = DefaultHasher::new();
    t.hash(&mut hasher);
    hasher.finish() as usize
}

/// Provides an initial label for an element based on a simple mapping (similar to using atomic numbers).
fn initial_label(element: &Element) -> usize {
    match element.kind {
        ElementType::H => 1,
        ElementType::C => 6,
        ElementType::N => 7,
        ElementType::O => 8,
        ElementType::F => 9,
        ElementType::P => 15,
        ElementType::S => 16,
        ElementType::Cl => 17,
        ElementType::Br => 35,
        ElementType::I => 53,
        ElementType::As => 33,
        ElementType::RGroup(n) => 1000 + n, // Arbitrary mapping for R-groups.
    }
}

/// Implements the Morgan algorithm for a given molecule graph.
/// 
/// # Arguments
/// * `graph` - A reference to the molecule graph (UnGraph<Element, Bond>).
/// * `max_iterations` - Maximum number of iterations to perform.
/// 
/// # Returns
/// A mapping from each node index to its final canonical label.
pub fn morgan_algorithm(graph: &UnGraph<Element, Bond>, max_iterations: usize) -> HashMap<NodeIndex, usize> {
    let mut labels: HashMap<NodeIndex, usize> = HashMap::new();
    
    // Initialize labels for all nodes.
    for node_ref in graph.node_references() {
        let node_index = node_ref.0;
        let element = node_ref.1;
        labels.insert(node_index, initial_label(element));
    }
    
    // Iteratively refine the labels.
    for _ in 0..max_iterations {
        let mut updated_labels = labels.clone();
        let mut changed = false;
        
        // Update each node's label based on its own label and its neighbors' labels.
        for node in graph.node_indices() {
            // Gather and sort the labels of neighboring nodes.
            let mut neighbor_labels: Vec<usize> = graph.neighbors(node)
                .map(|nbr| labels[&nbr])
                .collect();
            neighbor_labels.sort();
            
            // Combine current label with neighbor labels.
            let combined = (labels[&node], neighbor_labels);
            let new_label = compute_hash(&combined);
            
            // If the new label is different, record the change.
            if new_label != labels[&node] {
                updated_labels.insert(node, new_label);
                changed = true;
            }
        }
        
        labels = updated_labels;
        if !changed {
            break;
        }
    }
    
    labels
}

// Assume these types are imported from your module:
// use your_module::{Element, Bond, MoleculeGraph};

/// Rebuilds a canonical molecule graph from the Morgan algorithm labels.
/// 
/// # Arguments
/// * `graph` - A reference to the original molecule graph.
/// * `labels` - A mapping from each node's index to its canonical label.
/// 
/// # Returns
/// A new molecule graph with nodes ordered canonically.
pub fn rebuild_canonical_graph(
    graph: &UnGraph<Element, Bond>,
    labels: &HashMap<NodeIndex, usize>,
) -> UnGraph<Element, Bond> {
    // Create a vector of (old_node, label) pairs.
    let mut nodes_with_labels: Vec<(NodeIndex, usize)> = graph.node_indices()
        .map(|node| (node, *labels.get(&node).unwrap_or(&0)))
        .collect();

    // Sort nodes by canonical label, and use the node index as a tiebreaker.
    nodes_with_labels.sort_by_key(|(node, label)| (*label, node.index()));

    // Build a mapping from the old node index to a new node index.
    let mut mapping: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    let mut new_graph = UnGraph::<Element, Bond>::default();

    for (old_node, _) in nodes_with_labels.iter() {
        let new_index = new_graph.add_node(graph[*old_node]);
        mapping.insert(*old_node, new_index);
    }

    // Rebuild the edges using the new node indices.
    for edge in graph.edge_references() {
        let old_source = edge.source();
        let old_target = edge.target();
        let new_source = mapping[&old_source];
        let new_target = mapping[&old_target];
        new_graph.add_edge(new_source, new_target, *edge.weight());
    }

    new_graph
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn test_canonize() {
        let mol1 = parse_smiles("c1c(R)cccc1").unwrap();
        let mol2 = parse_smiles("c1ccccc1R").unwrap();
        let mol3 = parse_smiles("c1cccc(R)c1").unwrap();

        let canonized1 = mol1.morgan_canonize();
        let canonized2 = mol2.morgan_canonize();
        let canonized3 = mol3.morgan_canonize();
        
        let org1 = OrganicMolecule::from(mol1);
        org1.visualize("uncanonized1.png").unwrap();

        let org2 = OrganicMolecule::from(mol2);
        org2.visualize("uncanonized2.png").unwrap();

        let org3 = OrganicMolecule::from(mol3);
        org3.visualize("uncanonized3.png").unwrap();
        // let canonized1 = mol1.into_canon();
        // let canonized2 = mol2.into_canon();
        // let canonized1 = mol1.clone();
        // let canonized2 = mol2.clone();
        println!("Canonized 1: {:?}", canonized1);
        println!("Canonized 2: {:?}", canonized2);
        println!("Canonized 3: {:?}", canonized3);

        // println!("Is identical: {}", canonized1.is_identical(&canonized2));

        let org1 = OrganicMolecule::from(canonized1);
        org1.visualize("canonized1.png").unwrap();

        let org2 = OrganicMolecule::from(canonized2);
        org2.visualize("canonized2.png").unwrap();

        let org3 = OrganicMolecule::from(canonized3);
        org3.visualize("canonized3.png").unwrap();

        // Convert to smiles
        let smiles1 = org1.to_smiles().unwrap();
        let smiles2 = org2.to_smiles().unwrap();
        let smiles3 = org3.to_smiles().unwrap();

        println!("Smiles 1: {}", smiles1);
        println!("Smiles 2: {}", smiles2);
        println!("Smiles 3: {}", smiles3);

        assert_eq!(smiles1, smiles2);
        assert_eq!(smiles2, smiles3);
    }


    #[test]
    fn test_canonize2() {
        let mol1 = parse_smiles("c1c(R)cc(R)cc1").unwrap();
        let mol2 = parse_smiles("c1ccc(R)cc1R").unwrap();
        let mol3 = parse_smiles("c1cc(R)cc(R)c1").unwrap();
        
        // let canonized1 = mol1.into_canon();
        // let canonized2 = mol2.into_canon();
        // let canonized1 = mol1.clone();
        // let canonized2 = mol2.clone();
        let canonized1 = mol1.morgan_canonize();
        let canonized2 = mol2.morgan_canonize();
        let canonized3 = mol3.morgan_canonize();
        
        println!("Canonized 1: {:?}", canonized1);
        println!("Canonized 2: {:?}", canonized2);
        println!("Canonized 3: {:?}", canonized3);

        // println!("Is identical: {}", canonized1.is_identical(&canonized2));

        let org1 = OrganicMolecule::from(canonized1);
        org1.visualize("canonized1.png").unwrap();

        let org2 = OrganicMolecule::from(canonized2);
        org2.visualize("canonized2.png").unwrap();

        let org3 = OrganicMolecule::from(canonized3);
        org3.visualize("canonized2.png").unwrap();

        // Convert to smiles
        let smiles1 = org1.to_smiles().unwrap();
        let smiles2 = org2.to_smiles().unwrap();
        let smiles3 = org3.to_smiles().unwrap();

        println!("Smiles 1: {}", smiles1);
        println!("Smiles 2: {}", smiles2);
        println!("Smiles 3: {}", smiles3);

        assert_eq!(smiles1, smiles2);
        assert_eq!(smiles2, smiles3);
    }

    #[test]
    fn test_morgan_labels_consistency() {
        // Construct two molecule graphs representing the same molecule.
        let mol1 = parse_smiles("C1CC1").unwrap();
        let mol2 = parse_smiles("C1CC1").unwrap();
        let labels1 = morgan_algorithm(&mol1, 10);
        let labels2 = morgan_algorithm(&mol2, 10);
        
        // Since node indices may differ, compare sorted multisets of labels.
        let mut sorted_labels1: Vec<_> = labels1.values().cloned().collect();
        let mut sorted_labels2: Vec<_> = labels2.values().cloned().collect();
        sorted_labels1.sort();
        sorted_labels2.sort();

        assert!(is_identical(&mol1, &mol2), "Molecules should be identical");
        
        assert_eq!(sorted_labels1, sorted_labels2, "Morgan labels should be consistent for isomorphic graphs");
    }
    
    #[test]
    fn test_rebuild_canonical_graph_consistency() {
        // Test that rebuilding the graph using canonical labels preserves the SMILES representation.
        let mol = parse_smiles("c1ccccc1R").unwrap();
        let labels = morgan_algorithm(&mol, 10);
        let rebuilt = rebuild_canonical_graph(&mol, &labels);
        
        // Obtain the canonical graph via morgan_canonize, then compare SMILES from both representations.
        let canonized = mol.morgan_canonize();
        let org_original = OrganicMolecule::from(canonized);
        let org_rebuilt = OrganicMolecule::from(rebuilt);
        
        let smiles_original = org_original.to_smiles().unwrap();
        let smiles_rebuilt = org_rebuilt.to_smiles().unwrap();
        
        assert_eq!(smiles_original, smiles_rebuilt, "Rebuilt canonical graph should yield the same SMILES");
    }
    
    #[test]
    fn test_rebuild_graph_structure() {
        // Ensure that the rebuilt graph has the same number of nodes and edges as the original.
        let mol = parse_smiles("CCO").unwrap();
        let labels = morgan_algorithm(&mol, 10);
        let rebuilt = rebuild_canonical_graph(&mol, &labels);
        
        assert_eq!(mol.node_count(), rebuilt.node_count(), "Node count should be preserved");
        assert_eq!(mol.edge_count(), rebuilt.edge_count(), "Edge count should be preserved");
    }
    
    #[test]
    fn test_non_isomorphic_graphs() {
        // Test that non-isomorphic molecules (here, differing by a heteroatom) yield different Morgan labels.
        let mol1 = parse_smiles("CCO").unwrap();  // Ethanol
        let mol2 = parse_smiles("CCN").unwrap();  // Ethylamine, different heteroatom
        let labels1 = morgan_algorithm(&mol1, 10);
        let labels2 = morgan_algorithm(&mol2, 10);
        
        let mut sorted_labels1: Vec<_> = labels1.values().cloned().collect();
        let mut sorted_labels2: Vec<_> = labels2.values().cloned().collect();
        sorted_labels1.sort();
        sorted_labels2.sort();
        
        assert_ne!(sorted_labels1, sorted_labels2, "Non-isomorphic molecules should have different Morgan labels");
    }
}

/*
use petgraph::graph::{Graph, NodeIndex};
use petgraph::EdgeType;
use petgraph::visit::{EdgeRef, IntoNodeIdentifiers};
use std::collections::HashMap;

/// A trait for computing a canonical version of a structure.
/// The canonical form is an isomorphism invariant representation.
pub trait Canonize {
    /// Returns a canonical version of self.
    fn canonize(&self) -> Self;
}

/// Internal partition type used for refinement.
#[derive(Debug, Clone)]
struct Partition {
    blocks: Vec<Vec<NodeIndex>>,
}

impl Partition {
    fn new(blocks: Vec<Vec<NodeIndex>>) -> Self {
        Self { blocks }
    }

    /// Refine each block: for each node, build a signature (counts of neighbors in every block)
    /// and split the block by signature.
    fn refine<N, E, T>(&self, graph: &Graph<N, E, T>) -> Self
    where
        T: EdgeType,
    {
        let mut new_blocks = Vec::new();
        for block in &self.blocks {
            if block.len() <= 1 {
                new_blocks.push(block.clone());
            } else {
                let mut signature_map: HashMap<Vec<usize>, Vec<NodeIndex>> = HashMap::new();
                for &node in block {
                    let mut sig = Vec::new();
                    for other_block in &self.blocks {
                        let count = other_block.iter().filter(|&&other_node| {
                            graph.find_edge(node, other_node).is_some()
                                || graph.find_edge(other_node, node).is_some()
                        }).count();
                        sig.push(count);
                    }
                    signature_map.entry(sig).or_default().push(node);
                }
                for (_, nodes) in signature_map {
                    new_blocks.push(nodes);
                }
            }
        }
        Partition { blocks: new_blocks }
    }

    /// Iteratively refine until no further change occurs.
    fn stabilize<N, E, T>(&self, graph: &Graph<N, E, T>) -> Self
    where
        T: EdgeType,
    {
        let mut current = self.clone();
        while current.refine(graph).blocks != current.blocks {
            current = current.refine(graph);
        }
        current
    }
}

/// Returns true if every block in the partition is a singleton.
fn is_discrete(partition: &Partition) -> bool {
    partition.blocks.iter().all(|block| block.len() == 1)
}

/// Given a discrete (singleton) partition, produce a labeling vector mapping the
/// original node index to a canonical label (0..n-1). Sorting within blocks is done
/// by the node’s data (with a tie-breaker on the original index).
fn labeling_from_partition<N, E, T>(graph: &Graph<N, E, T>, partition: &Partition) -> Vec<usize>
where
    N: Ord + Clone,
    T: EdgeType,
{
    let mut blocks = partition.blocks.clone();
    for block in blocks.iter_mut() {
        block.sort_by_key(|node| (graph[*node].clone(), node.index()));
    }
    blocks.sort_by_key(|block| (graph[block[0]].clone(), block[0].index()));
    let flat: Vec<NodeIndex> = blocks.into_iter().flatten().collect();
    let mut labeling = vec![0; graph.node_count()];
    for (i, &node) in flat.iter().enumerate() {
        labeling[node.index()] = i;
    }
    labeling
}

/// Recursively compute a canonical labeling using individualization when needed.
fn compute_labeling<N, E, T>(
    graph: &Graph<N, E, T>,
    partition: Partition,
) -> Vec<usize>
where
    N: Ord + Clone,
    E: Ord + Clone,
    T: EdgeType,
{
    let refined = partition.stabilize(graph);
    if is_discrete(&refined) {
        return labeling_from_partition(graph, &refined);
    } else {
        // Find an ambiguous block (one with size > 1)
        let ambiguous_index = refined
            .blocks
            .iter()
            .position(|block| block.len() > 1)
            .unwrap();
        let ambiguous_block = refined.blocks[ambiguous_index].clone();

        let mut candidates: Vec<Vec<usize>> = Vec::new();
        // For each node in the ambiguous block, individualize it: force it into its own block.
        for &node in &ambiguous_block {
            let mut new_blocks = refined.blocks.clone();
            // Remove the ambiguous block.
            new_blocks.remove(ambiguous_index);
            // Add a singleton block for the chosen node.
            new_blocks.push(vec![node]);
            // And add a block with the remaining nodes (if any).
            let other: Vec<NodeIndex> = ambiguous_block.clone().into_iter().filter(|&x| x != node).collect();
            if !other.is_empty() {
                new_blocks.push(other);
            }
            let new_partition = Partition::new(new_blocks);
            let candidate = compute_labeling(graph, new_partition);
            candidates.push(candidate);
        }
        // Choose the lexicographically smallest labeling.
        candidates.sort();
        candidates[0].clone()
    }
}

/// Compute the canonical labeling of a graph.
fn canonical_labeling<N, E, T>(graph: &Graph<N, E, T>) -> Vec<usize>
where
    N: Ord + Clone,
    E: Ord + Clone,
    T: EdgeType,
{
    let mut groups: HashMap<usize, Vec<NodeIndex>> = HashMap::new();
    for node in graph.node_identifiers() {
        let degree = graph.edges(node).count();
        groups.entry(degree).or_default().push(node);
    }
    let initial = Partition::new(groups.into_iter().map(|(_, nodes)| nodes).collect());
    compute_labeling(graph, initial)
}

/// Implement the `Canonize` trait for petgraph’s `Graph`.
/// This routine uses the canonical labeling to rebuild the graph with nodes
/// labeled 0..n-1 and edges reinserted in a canonical (sorted) order.
impl<N, E, T> Canonize for Graph<N, E, T>
where
    N: Clone + Ord,
    E: Clone + Ord,
    T: EdgeType,
{
    fn canonize(&self) -> Self {
        let labeling = canonical_labeling(self);
        let n = self.node_count();

        // Build the inverse mapping: new label → original node index.
        let mut inverse_mapping = vec![0; n];
        for (orig, &new_label) in labeling.iter().enumerate() {
            inverse_mapping[new_label] = orig;
        }

        // Create a new graph by adding nodes in canonical order.
        let mut new_graph = Graph::<N, E, T>::with_capacity(n, self.edge_count());
        let mut new_node_indices = vec![NodeIndex::end(); n];
        for new_label in 0..n {
            let orig_index = inverse_mapping[new_label];
            let orig_node = NodeIndex::new(orig_index);
            let new_node = new_graph.add_node(self[orig_node].clone());
            new_node_indices[new_label] = new_node;
        }

        // Rebuild edges with canonical endpoints.
        let mut edges: Vec<(usize, usize, E)> = Vec::with_capacity(self.edge_count());
        for edge in self.edge_references() {
            let orig_source = edge.source();
            let orig_target = edge.target();
            let new_source = labeling[orig_source.index()];
            let new_target = labeling[orig_target.index()];
            let (s, t) = if self.is_directed() {
                (new_source, new_target)
            } else if new_source <= new_target {
                (new_source, new_target)
            } else {
                (new_target, new_source)
            };
            edges.push((s, t, edge.weight().clone()));
        }
        edges.sort_by(|(s1, t1, w1), (s2, t2, w2)| {
            s1.cmp(s2).then(t1.cmp(t2)).then(w1.cmp(w2))
        });
        for (s, t, weight) in edges {
            let source_node = new_node_indices[s];
            let target_node = new_node_indices[t];
            new_graph.add_edge(source_node, target_node, weight);
        }
        new_graph
    }
}


pub fn is_identical<N, E, T>(g1: &Graph<N, E, T>, g2: &Graph<N, E, T>) -> bool
where
    N: Ord + std::fmt::Debug,
    E: Ord + std::fmt::Debug,
    T: EdgeType,
{
    // Get all the node and edge indices
    let g1_nodes: Vec<_> = g1.node_indices().map(|n| &g2[n]).collect();
    let g2_nodes: Vec<_> = g2.node_indices().map(|n| &g2[n]).collect();

    let g1_edges: Vec<_> = g1.edge_references().map(|e| (e.source(), e.target())).collect();
    let g2_edges: Vec<_> = g2.edge_references().map(|e| (e.source(), e.target())).collect();

    // Check if the number of nodes and edges are the same
    if g1_nodes.len() != g2_nodes.len() || g1_edges.len() != g2_edges.len() {
        return false;
    }

    // Confirm all nodes are identical
    if g1_nodes != g2_nodes {
        return false;
    }

    if g1_edges != g2_edges {
        return false;
    }

    // Check if the node data is the same
    for (n1, n2) in g1_nodes.iter().zip(g2_nodes.iter()) {
        if n1 != n2 {
            return false;
        }
    }

    // Check if the edge data is the same
    for (e1, e2) in g1_edges.iter().zip(g2_edges.iter()) {
        if e1 != e2 {
            return false;
        }
    }

    true
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;
    use petgraph::Graph;

    #[test]
    fn test_canonize_cycle() {
        // Create a 4-node cycle graph.
        let mut graph: Graph<&str, &str, _> = Graph::new_undirected();
        let a = graph.add_node("a");
        let b = graph.add_node("b");
        let c = graph.add_node("c");
        let d = graph.add_node("d");
        graph.add_edge(a, b, "ab");
        graph.add_edge(b, c, "bc");
        graph.add_edge(c, d, "cd");
        graph.add_edge(d, a, "da");

        let canon = graph.canonize();
        // The canonical graph should have 4 nodes (with new indices 0,1,2,3)
        // and the edges reinserted in a fixed order.
        assert_eq!(canon.node_count(), 4);
        assert_eq!(canon.edge_count(), 4);
    }

    #[test]
    fn test_isomorphic_canonization() {
        // Construct two isomorphic 4-cycle graphs with different node insertion orders.
        let mut graph1: Graph<&str, &str, _> = Graph::new_undirected();
        let a1 = graph1.add_node("a");
        let b1 = graph1.add_node("b");
        let c1 = graph1.add_node("c");
        let d1 = graph1.add_node("d");
        graph1.add_edge(a1, b1, "ab");
        graph1.add_edge(b1, c1, "bc");
        graph1.add_edge(c1, d1, "cd");
        graph1.add_edge(d1, a1, "da");

        let mut graph2: Graph<&str, &str, _> = Graph::new_undirected();
        let d2 = graph2.add_node("d");
        let c2 = graph2.add_node("c");
        let b2 = graph2.add_node("b");
        let a2 = graph2.add_node("a");
        graph2.add_edge(a2, b2, "ab");
        graph2.add_edge(b2, c2, "bc");
        graph2.add_edge(c2, d2, "cd");
        graph2.add_edge(d2, a2, "da");

        let canon1 = graph1.canonize();
        let canon2 = graph2.canonize();
        
        // The canonical forms should be identical.
        assert!(is_identical(&canon1, &canon2));
        // The original graphs should not be identical.
        assert!(!is_identical(&graph1, &graph2));
    }

    #[test]
    fn test_isomorphic_canonization2() {
        // Construct two isomorphic 4-cycle graphs with different node insertion orders.
        let mut graph1: Graph<&str, &str, _> = Graph::new_undirected();
        let a1 = graph1.add_node("a");
        let b1 = graph1.add_node("b");
        let c1 = graph1.add_node("c");
        let d1 = graph1.add_node("d");
        let e1 = graph1.add_node("e");
        graph1.add_edge(a1, b1, "ab");
        graph1.add_edge(b1, c1, "bc");
        graph1.add_edge(c1, d1, "cd");
        graph1.add_edge(d1, a1, "da");
        graph1.add_edge(a1, e1, "ae");

        let mut graph2: Graph<&str, &str, _> = Graph::new_undirected();
        let e2 = graph2.add_node("e");
        let d2 = graph2.add_node("d");
        let c2 = graph2.add_node("c");
        let b2 = graph2.add_node("b");
        let a2 = graph2.add_node("a");
        graph2.add_edge(a2, b2, "ab");
        graph2.add_edge(b2, c2, "bc");
        graph2.add_edge(c2, d2, "cd");
        graph2.add_edge(d2, a2, "da");
        graph2.add_edge(a2, e2, "ae");


        let canon1 = graph1.canonize();
        let canon2 = graph2.canonize();
        println!("Graph 1: {:?}", graph1);
        println!("Graph 2: {:?}", graph2);

        println!("Canon 1: {:?}", canon1);
        println!("Canon 2: {:?}", canon2);
        
        // The canonical forms should be identical.
        assert!(is_identical(&canon1, &canon2));
        // The original graphs should not be identical.
        assert!(!is_identical(&graph1, &graph2));
    }

    /// Test that two isomorphic graphs produce the same canonical fingerprint.
    #[test]
    fn test_isomorphic_graphs() {
        // Construct first graph: a 4-node cycle.
        let mut graph1: Graph<(), (), _> = Graph::new_undirected();
        let a1 = graph1.add_node(()); // will be 0
        let b1 = graph1.add_node(()); // 1
        let c1 = graph1.add_node(()); // 2
        let d1 = graph1.add_node(()); // 3

        graph1.add_edge(a1, b1, ());
        graph1.add_edge(b1, c1, ());
        graph1.add_edge(c1, d1, ());
        graph1.add_edge(d1, a1, ());

        // Construct second graph: also a 4-node cycle but with a different insertion order.
        let mut graph2: Graph<(), (), _> = Graph::new_undirected();
        let d2 = graph2.add_node(()); // inserted first
        let c2 = graph2.add_node(());
        let b2 = graph2.add_node(());
        let a2 = graph2.add_node(()); // inserted last

        graph2.add_edge(a2, b2, ());
        graph2.add_edge(b2, c2, ());
        graph2.add_edge(c2, d2, ());
        graph2.add_edge(d2, a2, ());
        

        let canon1 = graph1.canonize();
        let canon2 = graph2.canonize();
        println!("Graph 1: {:?}", canon1);
        println!("Graph 2: {:?}", graph2);

        println!("Canon 1: {:?}", canon1);
        println!("Canon 2: {:?}", canon2);
        assert!(!is_identical(&graph1, &graph2));
        assert!(is_identical(&canon1, &canon2));
    }

    /// Test that two isomorphic graphs produce the same canonical fingerprint.
    #[test]
    fn test_isomorphic_graphs2() {
        // Construct first graph: a 4-node cycle.
        let mut graph1: Graph<(), (), _> = Graph::new_undirected();
        let a1 = graph1.add_node(()); // will be 0
        let b1 = graph1.add_node(()); // 1
        let c1 = graph1.add_node(()); // 2
        let d1 = graph1.add_node(()); // 3
        let e1 = graph1.add_node(()); // 4

        graph1.add_edge(a1, b1, ());
        graph1.add_edge(b1, c1, ());
        graph1.add_edge(c1, d1, ());
        graph1.add_edge(d1, a1, ());
        graph1.add_edge(d1, e1, ());

        // Construct second graph: also a 4-node cycle but with a different insertion order.
        let mut graph2: Graph<(), (), _> = Graph::new_undirected();
        let e2 = graph2.add_node(()); // inserted first
        let d2 = graph2.add_node(()); // inserted first
        let c2 = graph2.add_node(());
        let b2 = graph2.add_node(());
        let a2 = graph2.add_node(()); // inserted last

        graph2.add_edge(a2, b2, ());
        graph2.add_edge(b2, c2, ());
        graph2.add_edge(c2, d2, ());
        graph2.add_edge(d2, a2, ());
        graph2.add_edge(a2, b2, ());
        

        let canon1 = graph1.canonize();
        let canon2 = graph2.canonize();
        println!("Graph 1: {:?}", canon1);
        println!("Graph 2: {:?}", graph2);

        println!("Canon 1: {:?}", canon1);
        println!("Canon 2: {:?}", canon2);
        assert!(!is_identical(&graph1, &graph2));
        assert!(is_identical(&canon1, &canon2));
    }

    #[test]
    fn test_canonize_smiles() {
        let smiles1 = "c1cc(R)ccc1"; // Benzene
        let smiles2 = "c1ccccc1R"; // Benzene
        let mol1 = canonize_smiles(smiles1).unwrap();
        let mol2 = canonize_smiles(smiles2).unwrap();

        // let graph = molecule_graph_to_unit_graph(&mol);
        // let canon = graph.canonize();
        println!("Canonical SMILES #1: {:?}", mol1);
        println!("Canonical SMILES #2: {:?}", mol2);

        assert_ne!(smiles1, smiles2);
        assert_eq!(mol1, mol2);
    }
}
*/

/*
use petgraph::graph::{Graph, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use petgraph::{EdgeType, stable_graph::IndexType};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use tracing::{debug, trace, instrument};

#[derive(Debug, Clone)]
pub struct OrbitSet {
    parent: Vec<usize>,
}

impl OrbitSet {
    fn new(n: usize) -> Self {
        OrbitSet {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] == x {
            x
        } else {
            let p = self.parent[x];
            let r = self.find(p);
            self.parent[x] = r;
            r
        }
    }
    fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx != ry {
            self.parent[ry] = rx;
        }
    }
    pub fn same_orbit(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }
}

#[derive(Debug, Clone)]
pub struct SearchState<Ix> {
    base_labeling: Vec<NodeIndex<Ix>>,
    pos_of_node: Vec<usize>,
    pub orbits: OrbitSet,

    /// The partial labeling array we reorder. 
    /// If `depth = k`, then `working_labeling[0..k]` is fixed, 
    /// and `[k..]` might be unassigned or temporarily swapped in.
    pub working_labeling: Vec<NodeIndex<Ix>>,
    pub depth: usize,

    /// The best final labeling so far (global minimum).
    pub best_labeling: Option<Vec<NodeIndex<Ix>>>,
    
    /// For partial check: partition states we've seen at each depth -> record a "representative" partial labeling
    pub seen_partitions: HashMap<(usize, PartitionSignature), Vec<NodeIndex<Ix>>>,
}

/// A small signature capturing the partition “shape” for partial checks:
/// e.g. a sorted list of cell sizes, or something that identifies adjacency splits, etc.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PartitionSignature {
    /// We just store each cell’s sorted node indices. 
    /// For partial checks, you might want something smaller, but let's keep it simple.
    cells: Vec<Vec<usize>>,
}

/// Our main entry
pub fn canonical_labeling<N, E, Ty, Ix>(graph: &Graph<N, E, Ty, Ix>) -> Vec<NodeIndex<Ix>>
where
    N: Ord + Clone,
    E: Clone,
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut base_label: Vec<NodeIndex<Ix>> = graph.node_indices().collect();
    base_label.sort_by_key(|n| n.index());
    let n = base_label.len();

    // Build initial partition
    let mut partition = initial_partition_label_degree(graph, &base_label);
    while refine_equitable(graph, &mut partition) {}

    let mut st = SearchState {
        base_labeling: base_label.clone(),
        pos_of_node: vec![0; n],
        orbits: OrbitSet::new(n),
        working_labeling: base_label.clone(),
        depth: 0,
        best_labeling: None,
        seen_partitions: HashMap::new(),
    };
    // fill pos_of_node
    for (i, &nd) in st.base_labeling.iter().enumerate() {
        st.pos_of_node[nd.index()] = i;
    }

    backtrack_refine(graph, &mut partition, &mut st);
    st.best_labeling.unwrap_or(st.working_labeling)
}

// -----------------------------------------------------------------------------
// Partition initialization & Equitable refinement
// -----------------------------------------------------------------------------

fn initial_partition_label_degree<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    base_labeling: &[NodeIndex<Ix>],
) -> Vec<Vec<NodeIndex<Ix>>>
where
    N: Ord + Clone,
    Ty: EdgeType,
    Ix: IndexType,
{
    use std::collections::BTreeMap;
    let mut map: BTreeMap<(N, usize), Vec<NodeIndex<Ix>>> = BTreeMap::new();
    for &v in base_labeling {
        let lbl = graph[v].clone();
        let deg = graph.edges(v).count();
        map.entry((lbl, deg)).or_default().push(v);
    }
    let mut partition: Vec<Vec<NodeIndex<Ix>>> = map.into_values().collect();
    for c in &mut partition {
        c.sort_by_key(|n| n.index());
    }
    partition.sort_by_key(|c| c[0].index());
    partition
}

fn refine_equitable<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    partition: &mut Vec<Vec<NodeIndex<Ix>>>,
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut changed = false;
    let mut new_partition = Vec::new();

    for cell in partition.drain(..) {
        if cell.len() <= 1 {
            new_partition.push(cell);
            continue;
        }
        use std::collections::BTreeMap;
        let mut groups: BTreeMap<Vec<usize>, Vec<NodeIndex<Ix>>> = BTreeMap::new();
        for &v in &cell {
            let mut sig = Vec::with_capacity(new_partition.len() + 1);
            for group in &new_partition {
                let count = group.iter().filter(|&&w| is_adj(graph, v, w)).count();
                sig.push(count);
            }
            let self_count = cell.iter().filter(|&&w| w != v && is_adj(graph, v, w)).count();
            sig.push(self_count);

            groups.entry(sig).or_default().push(v);
        }
        if groups.len() > 1 {
            changed = true;
        }
        let mut subcells: Vec<Vec<NodeIndex<Ix>>> = groups.into_values().collect();
        for sc in &mut subcells {
            sc.sort_by_key(|n| n.index());
        }
        subcells.sort_by_key(|sc| sc[0].index());
        new_partition.extend(subcells);
    }
    new_partition.sort_by_key(|c| c[0].index());
    *partition = new_partition;
    changed
}

// -----------------------------------------------------------------------------
// Partial Backtracking with Intra-level partial checks
// -----------------------------------------------------------------------------

#[instrument(skip_all, fields(depth = st.depth, part_len = partition.len()))]
fn backtrack_refine<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    partition: &mut Vec<Vec<NodeIndex<Ix>>>,
    st: &mut SearchState<Ix>,
) where
    N: Ord + Clone,
    E: Clone,
    Ty: EdgeType,
    Ix: IndexType,
{
    let n = st.base_labeling.len();
    // If discrete => full labeling
    if partition.iter().all(|c| c.len() == 1) {
        // unify orbits if automorphism
        unify_if_automorphism(graph, &st.working_labeling.clone(), st);

        // check if better => update best_labeling
        match &st.best_labeling {
            None => st.best_labeling = Some(st.working_labeling.clone()),
            Some(cur) => {
                if is_labeling_better(graph, &st.working_labeling, cur) {
                    st.best_labeling = Some(st.working_labeling.clone());
                }
            }
        }
        return;
    }

    // Partial check: see if this partition shape has been encountered at st.depth.
    // If so, check if the partial labeling is an automorphism => unify orbits => maybe skip
    let part_sig = partition_signature(partition);
    if let Some(prev_lab) = st.seen_partitions.get(&(st.depth, part_sig.clone())).cloned() {
        // we have a "representative" partial labeling that produced this partition shape
        // check if the current partial labeling is an automorphism of that one
        if partial_is_automorphism(graph, &st.base_labeling, &prev_lab, &st.working_labeling, st.depth) {
            // unify orbits among positions in [0..depth], skip further expansions
            unify_partial_orbits(&prev_lab, &st.working_labeling.clone(), st.depth, st);
            // This may unify orbits so that we skip symmetrical pivots 
            // but for a complete skip, you'd do it if the partial labeling is *identical* or something 
            // We can do an early return if we consider them fully identical 
            // (though that might skip distinct expansions). 
            // Let's just unify orbits and *not* forcibly skip. 
            // Real code might do more advanced logic:
            return;
        }
    } else {
        // store this partial labeling as a representative
        st.seen_partitions.insert((st.depth, part_sig.clone()), st.working_labeling.clone());
    }

    // Not discrete => pick target cell, pivot each node not in same orbit
    let cell_idx = pick_target_cell(partition);
    let mut cell_nodes = partition[cell_idx].clone();
    cell_nodes.sort_by_key(|n| n.index());

    let old_depth = st.depth;
    for (i, &pivot_node) in cell_nodes.iter().enumerate() {
        // Orbit skip
        let pivot_pos = st.pos_of_node[pivot_node.index()];
        let mut skip = false;
        for &prev in &cell_nodes[..i] {
            let pp = st.pos_of_node[prev.index()];
            if st.orbits.same_orbit(pp, pivot_pos) {
                skip = true;
                break;
            }
        }
        if skip {
            continue;
        }

        // assign pivot_node to st.working_labeling[depth]
        let d = st.depth;
        let old_node = st.working_labeling[d];
        if old_node != pivot_node {
            let pivot_idx = st.working_labeling.iter().position(|&x| x == pivot_node).unwrap();
            st.working_labeling.swap(d, pivot_idx);
        }

        // remove pivot_node from partition cell
        let saved_partition = partition.clone();
        partition[cell_idx].retain(|&x| x != pivot_node);
        partition.push(vec![pivot_node]);
        partition.sort_by_key(|c| c[0].index());

        // refine
        while refine_equitable(graph, partition) {}

        // recursion deeper
        st.depth += 1;
        backtrack_refine(graph, partition, st);
        st.depth = d;

        // revert partition
        *partition = saved_partition;

        // revert st.working_labeling
        let pivot_now = st.working_labeling[d];
        if pivot_now != old_node {
            let idx_now = st.working_labeling.iter().position(|&x| x == pivot_now).unwrap();
            let old_idx = st.working_labeling.iter().position(|&x| x == old_node).unwrap();
            st.working_labeling.swap(idx_now, old_idx);
        }
    }
    st.depth = old_depth;
}

/// Build a “partition signature” for partial checks. 
/// For example, store each cell’s sorted node.index().
fn partition_signature<Ix: IndexType>(partition: &[Vec<NodeIndex<Ix>>]) -> PartitionSignature {
    let mut cells: Vec<Vec<usize>> = partition
        .iter()
        .map(|c| c.iter().map(|x| x.index()).collect())
        .collect();
    for c in &mut cells {
        c.sort();
    }
    cells.sort_by_key(|c| c[0]);
    PartitionSignature { cells }
}

/// If `candidate` is an automorphism of `base_lab` on the first `depth` positions, unify orbits.
fn unify_if_automorphism<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    candidate: &[NodeIndex<Ix>],
    st: &mut SearchState<Ix>,
) where
    N: Ord + Clone,
    E: Clone,
    Ty: EdgeType,
    Ix: IndexType,
{
    let n = st.base_labeling.len();
    if is_automorphism_of_base(graph, candidate, &st.base_labeling) {
        // unify all positions
        let mut rev_base = BTreeMap::new();
        for (j, &bn) in st.base_labeling.iter().enumerate() {
            rev_base.insert(bn, j);
        }
        for i in 0..n {
            let mapped_node = candidate[i];
            let j = rev_base[&mapped_node];
            st.orbits.union(i, j);
        }
        debug!("Found automorphism => unify orbits");
    }
}

/// Partial isomorphism check: Are the *first `depth` positions* in `candB` an automorphism of `candA`?
/// i.e. do they preserve adjacency among those assigned nodes?
fn partial_is_automorphism<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    base: &[NodeIndex<Ix>],
    candA: &[NodeIndex<Ix>],
    candB: &[NodeIndex<Ix>],
    depth: usize,
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    // If depth=0 or 1, trivial. If depth=2 or more, let's check adjacency among candA[0..depth].
    if depth < 2 { return false; } 
    // build adjacency sets for candB among first `depth`.
    let mut edges_b = BTreeSet::new();
    for i in 0..depth {
        for j in (i+1)..depth {
            if is_adj(graph, candB[i], candB[j]) {
                edges_b.insert((i, j));
            }
        }
    }
    // Check adjacency in A
    for i in 0..depth {
        for j in (i+1)..depth {
            let hasA = is_adj(graph, candA[i], candA[j]);
            let hasB = edges_b.contains(&(i, j));
            if hasA != hasB {
                return false;
            }
        }
    }
    true
}

/// If partial labeling candB is automorphism of candA up to `depth`, unify orbits among the positions in [0..depth].
fn unify_partial_orbits<Ix: IndexType>(
    candA: &[NodeIndex<Ix>],
    candB: &[NodeIndex<Ix>],
    depth: usize,
    st: &mut SearchState<Ix>,
) {
    // E.g. if candA[i] == nodeX, candB[i] == nodeY, then unify orbits of X and Y in union-find.
    // But we must unify the *positions* in base_labeling:
    //    posA = st.pos_of_node[nodeX.index()]
    //    posB = st.pos_of_node[nodeY.index()]
    for i in 0..depth {
        let x = candA[i];
        let y = candB[i];
        let px = st.pos_of_node[x.index()];
        let py = st.pos_of_node[y.index()];
        st.orbits.union(px, py);
    }
    debug!("Unified orbits from partial automorphism among first {} positions", depth);
}

// -----------------------------------------------------------------------------
// Automorphism checks, labeling comparison, etc.
// -----------------------------------------------------------------------------

fn is_automorphism_of_base<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    candidate: &[NodeIndex<Ix>],
    base_label: &[NodeIndex<Ix>],
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let n = base_label.len();
    if candidate.len() != n {
        return false;
    }
    let mut edgesetB = BTreeSet::new();
    for i in 0..n {
        for j in (i+1)..n {
            if is_adj(graph, candidate[i], candidate[j]) {
                edgesetB.insert((i, j));
            }
        }
    }
    for i in 0..n {
        for j in (i+1)..n {
            let has_base = is_adj(graph, base_label[i], base_label[j]);
            let has_cand = edgesetB.contains(&(i, j));
            if has_base != has_cand {
                return false;
            }
        }
    }
    true
}

fn is_labeling_better<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    new_lab: &[NodeIndex<Ix>],
    cur_lab: &[NodeIndex<Ix>],
) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let sig_new = candidate_signature(graph, new_lab);
    let sig_cur = candidate_signature(graph, cur_lab);
    match sig_new.cmp(&sig_cur) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => compare_labelings(new_lab, cur_lab) == Ordering::Less,
    }
}

fn candidate_signature<N, E, Ty, Ix>(
    graph: &Graph<N, E, Ty, Ix>,
    labeling: &[NodeIndex<Ix>],
) -> Vec<(usize, usize)>
where
    Ty: EdgeType,
    Ix: IndexType,
{
    let mut pos = BTreeMap::new();
    for (i, &nd) in labeling.iter().enumerate() {
        pos.insert(nd.index(), i);
    }
    let mut edges = Vec::new();
    for e in graph.edge_references() {
        let s = pos[&e.source().index()];
        let t = pos[&e.target().index()];
        edges.push((s.min(t), s.max(t)));
    }
    edges.sort();
    edges
}

fn compare_labelings<Ix: IndexType>(
    a: &[NodeIndex<Ix>],
    b: &[NodeIndex<Ix>],
) -> Ordering {
    for i in 0..a.len().min(b.len()) {
        let cmp = a[i].index().cmp(&b[i].index());
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    a.len().cmp(&b.len())
}

fn is_adj<N, E, Ty, Ix>(graph: &Graph<N, E, Ty, Ix>, a: NodeIndex<Ix>, b: NodeIndex<Ix>) -> bool
where
    Ty: EdgeType,
    Ix: IndexType,
{
    graph.find_edge(a, b).is_some() || graph.find_edge(b, a).is_some()
}

/// Return the index of the target cell
fn pick_target_cell<Ix: IndexType>(partition: &[Vec<NodeIndex<Ix>>]) -> usize {
    let mut chosen = None;
    let mut best_min = None;
    for (i, c) in partition.iter().enumerate() {
        if c.len() > 1 {
            let m = c.iter().map(|n| n.index()).min().unwrap();
            if best_min.is_none() || m < best_min.unwrap() {
                best_min = Some(m);
                chosen = Some(i);
            }
        }
    }
    chosen.expect("No non-singleton cell found, but partition not discrete?")
}

// -----------------------------------------------------------------------------
// "canonize" to build a new graph from the final labeling
// -----------------------------------------------------------------------------

pub fn canonize<N, E, Ix>(g: &UnGraph<N, E, Ix>) -> UnGraph<N, E, Ix>
where
    N: Ord + Clone,
    E: Ord + Clone,
    Ix: IndexType,
{
    let labeling = canonical_labeling(g);
    let mut new_graph = UnGraph::default();
    let mut map: BTreeMap<NodeIndex<Ix>, NodeIndex<Ix>> = BTreeMap::new();
    for &old in &labeling {
        let data = g[old].clone();
        let new_idx = new_graph.add_node(data);
        map.insert(old, new_idx);
    }
    for e in g.edge_references() {
        let s = map[&e.source()];
        let t = map[&e.target()];
        new_graph.add_edge(s, t, e.weight().clone());
    }
    new_graph
}


#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::Graph;
    use crate::*;

    #[test]
    fn test_small_example() {
        let mut g = UnGraph::new_undirected();
        let a = g.add_node("a");
        let b = g.add_node("b");
        let c = g.add_node("c");
        let d = g.add_node("d");

        g.add_edge(a, b, "");
        g.add_edge(b, c, "");
        g.add_edge(c, d, "");
        g.add_edge(d, a, "");

        let labeling = canonical_labeling(&g);
        // Check it's a permutation of [a,b,c,d].
        assert_eq!(labeling.len(), 4);

        // For manual debugging, print it out:
        println!("Labeling for small example: {:?}", labeling);
    }

    #[test]
    fn test_isomorphic_cycles() {
        init_logging("trace");
        // Two isomorphic cycles, different node insertion order.
        let mut g1 = UnGraph::new_undirected();
        let a1 = g1.add_node("a1");
        let b1 = g1.add_node("b1");
        let c1 = g1.add_node("c1");
        let d1 = g1.add_node("d1");
        g1.add_edge(a1, b1, "");
        g1.add_edge(b1, c1, "");
        g1.add_edge(c1, d1, "");
        g1.add_edge(d1, a1, "");

        let mut g2 = UnGraph::new_undirected();
        let d2 = g2.add_node("d2");
        let a2 = g2.add_node("a2");
        let c2 = g2.add_node("c2");
        let b2 = g2.add_node("b2");
        g2.add_edge(d2, a2, "");
        g2.add_edge(a2, b2, "");
        g2.add_edge(b2, c2, "");
        g2.add_edge(c2, d2, "");

        let lab1 = canonical_labeling(&g1);
        let lab2 = canonical_labeling(&g2);

        // Convert each labeling to a canonical edge form to compare.
        let cf1 = canonical_edge_form(&g1, &lab1);
        let cf2 = canonical_edge_form(&g2, &lab2);

        assert_eq!(
            cf1, cf2,
            "Isomorphic cycles should yield the same canonical edge form"
        );
    }

    // Helper to convert a labeling into a sorted edge form.
    fn canonical_edge_form<N, E, Ty, Ix>(
        graph: &Graph<N, E, Ty, Ix>,
        labeling: &[NodeIndex<Ix>],
    ) -> Vec<(usize, usize)>
    where
        Ty: EdgeType,
        Ix: petgraph::graph::IndexType,
    {
        // Build mapping from old index -> new position in labeling.
        let mut map = HashMap::new();
        for (new_pos, &old_idx) in labeling.iter().enumerate() {
            map.insert(old_idx, new_pos);
        }

        let mut edges = Vec::new();
        for edge in graph.edge_references() {
            let s_new = map[&edge.source()];
            let t_new = map[&edge.target()];
            let (u, v) = if graph.is_directed() {
                (s_new, t_new)
            } else {
                (s_new.min(t_new), s_new.max(t_new))
            };
            edges.push((u, v));
        }
        edges.sort();
        edges
    }

    #[test]
    fn test_canonize() {
        use nauty_pet::prelude::*;
        use nauty_pet::canon::*;
        let mol1 = parse_smiles("c1c(R)cccc1").unwrap();
        let mol2 = parse_smiles("c1ccccc1R").unwrap();
        let mol3 = parse_smiles("c1cccc(R)c1").unwrap();
        
        // let canonized1 = mol1.into_canon();
        // let canonized2 = mol2.into_canon();
        // let canonized1 = mol1.clone();
        // let canonized2 = mol2.clone();
        let canonized1 = canonize(&mol1);
        let canonized2 = canonize(&mol2);
        let canonized3 = canonize(&mol3);
        
        println!("Canonized 1: {:?}", canonized1);
        println!("Canonized 2: {:?}", canonized2);
        println!("Canonized 3: {:?}", canonized3);

        // println!("Is identical: {}", canonized1.is_identical(&canonized2));

        let org1 = OrganicMolecule::from(canonized1);
        org1.visualize("canonized1.png").unwrap();

        let org2 = OrganicMolecule::from(canonized2);
        org2.visualize("canonized2.png").unwrap();

        let org3 = OrganicMolecule::from(canonized3);
        org3.visualize("canonized2.png").unwrap();

        // Convert to smiles
        let smiles1 = org1.to_smiles().unwrap();
        let smiles2 = org2.to_smiles().unwrap();
        let smiles3 = org3.to_smiles().unwrap();

        println!("Smiles 1: {}", smiles1);
        println!("Smiles 2: {}", smiles2);
        println!("Smiles 3: {}", smiles3);
    }
}


#[cfg(test)]
mod unit_tests {
    use super::*;
    use crate::*;
    use petgraph::prelude::*;
    use petgraph::Graph;
    use std::collections::HashMap;

    #[test]
    fn test_candidate_signature() {
        // Let's build a small graph:
        //   0 -- 1
        //   |
        //   2
        // Then define a labeling and check the edge signature we get.

        let mut g = UnGraph::<(), ()>::new_undirected();
        let n0 = g.add_node(()); // index 0
        let n1 = g.add_node(()); // index 1
        let n2 = g.add_node(()); // index 2
        g.add_edge(n0, n1, ());
        g.add_edge(n0, n2, ());

        // Let’s say the labeling is [n0, n1, n2].
        // Then in the new order: node0 -> 0, node1 -> 1, node2 -> 2.
        // We have edges (0,1) and (0,2) in the sorted signature.
        let label = vec![n0, n1, n2];
        let sig = candidate_signature(&g, &label);

        assert_eq!(sig, vec![(0,1), (0,2)]);
    }

    #[test]
    fn test_compare_labelings() {
        // We'll build two labelings: [0,1,2] and [0,2,1] and compare them
        // The first difference is at index 1: 1 < 2, so [0,1,2] is "less".

        let lab_a: Vec<NodeIndex<u32>> = vec![NodeIndex::new(0), NodeIndex::new(1), NodeIndex::new(2)];
        let lab_b = vec![NodeIndex::new(0), NodeIndex::new(2), NodeIndex::new(1)];

        let ordering = compare_labelings(&lab_a, &lab_b);
        assert_eq!(ordering, std::cmp::Ordering::Less);

        let ordering_rev = compare_labelings(&lab_b, &lab_a);
        assert_eq!(ordering_rev, std::cmp::Ordering::Greater);

        // If they are equal, we return Ordering::Equal
        let ordering_same = compare_labelings(&lab_a, &lab_a);
        assert_eq!(ordering_same, std::cmp::Ordering::Equal);
    }

    #[test]
    fn test_canonical_labeling() {
        // For a simple path of length 3, we can see what the canonical labeling does.
        // All that matters is that the code picks one consistent labeling for the path.
        let mut g = UnGraph::<(), ()>::new_undirected();
        let n0 = g.add_node(()); // index 0
        let n1 = g.add_node(()); // index 1
        let n2 = g.add_node(()); // index 2
        g.add_edge(n0, n1, ());
        g.add_edge(n1, n2, ());

        let labeling = canonical_labeling(&g);

        // For a path, you might get [0,1,2] or [2,1,0], etc.
        // There's no single "right" answer, but we want to confirm it’s a valid permutation
        // and it's stable for repeated calls (i.e., the function is deterministic).
        let labeling2 = canonical_labeling(&g);

        assert_eq!(labeling.len(), 3);
        assert_eq!(labeling2, labeling, "Canonical labeling should be deterministic");
    }

    #[test]
    fn test_canonize() {
        // We'll use a cycle of 4. Then reorder the nodes in the input graph, call canonize,
        // and ensure we get the same final graph. This checks isomorphism invariance.
        init_logging("trace");
        let mut g1 = UnGraph::<char, ()>::new_undirected();
        let a1 = g1.add_node('a');
        let b1 = g1.add_node('b');
        let c1 = g1.add_node('c');
        let d1 = g1.add_node('d');
        g1.add_edge(a1, b1, ());
        g1.add_edge(b1, c1, ());
        g1.add_edge(c1, d1, ());
        g1.add_edge(d1, a1, ());

        // Now create the same cycle but with a different insertion order
        let mut g2 = UnGraph::<char, ()>::new_undirected();
        let d2 = g2.add_node('d');
        let a2 = g2.add_node('a');
        let c2 = g2.add_node('c');
        let b2 = g2.add_node('b');
        g2.add_edge(d2, a2, ());
        g2.add_edge(a2, b2, ());
        g2.add_edge(b2, c2, ());
        g2.add_edge(c2, d2, ());
        info!("Canonizing Graph 1: {:?}", g1);
        let c1 = canonize(&g1);
        info!("Canonizing Graph 2: {:?}", g2);
        let c2 = canonize(&g2);
        // For two isomorphic graphs, canonize should produce an identical “shape”.
        // This is a bit tricky to test for equality because they might have
        // different internal node indices but the same structure. We can do
        // an edge-by-edge check or a “graph isomorphism check.” For simplicity, compare:
        assert_eq!(c1.node_count(), c2.node_count());
        assert_eq!(c1.edge_count(), c2.edge_count());

        // Or do a quick adjacency check in a simple way:
        // There's an even simpler approach if you trust that `canonize` places the
        // final nodes in a consistent order, so let's do an edge signature check:
        let sig1 = candidate_signature(&c1, &c1.node_indices().collect::<Vec<_>>());
        let sig2 = candidate_signature(&c2, &c2.node_indices().collect::<Vec<_>>());
        assert_eq!(sig1, sig2, "Canonized results should have the same edge signature");
    }
}



#[cfg(test)]
mod harder_tests {
    use crate::init_logging;

    use super::*;
    use petgraph::prelude::*;
    use petgraph::Graph;

    /// A utility that, given a set of graphs, canonizes them all
    /// and asserts they produce the same canonical edge signature.
    /// Useful for testing multiple isomorphic constructions.
    fn assert_isomorphic_canon<N, E>(graphs: Vec<UnGraph<N, E>>)
    where
        N: Clone + Ord,
        E: Clone + Ord,
    {
        init_logging("trace");
        let canonized: Vec<UnGraph<N, E>> =
            graphs.into_iter().map(|g| canonize(&g.into())).collect();

        // Compare each canonized graph’s edge signature to the first one.
        let sig0 = candidate_signature(
            &canonized[0],
            &canonized[0].node_indices().collect::<Vec<_>>(),
        );

        for c in &canonized[1..] {
            let sig = candidate_signature(c, &c.node_indices().collect::<Vec<_>>());
            assert_eq!(
                sig0, sig,
                "All isomorphic graphs should have the same signature after canonization"
            );
        }
    }

    #[test]
    fn test_complete_graph_k5() {
        // K5: complete graph on 5 nodes. This is a highly symmetric graph.
        // We construct it in two different orders but it should canonize to the same form.

        // Graph 1
        let mut g1 = UnGraph::<usize, ()>::new_undirected();
        let nodes1 = (0..5).map(|i| g1.add_node(i)).collect::<Vec<_>>();
        // Add edges for all pairs
        for i in 0..5 {
            for j in (i + 1)..5 {
                g1.add_edge(nodes1[i], nodes1[j], ());
            }
        }

        // Graph 2, random insertion order
        let mut g2 = UnGraph::<usize, ()>::new_undirected();
        // Shuffle node creation
        let order = vec![2, 4, 1, 0, 3];
        let mut map2 = Vec::with_capacity(5);
        for &i in &order {
            map2.push((i, g2.add_node(i)));
        }
        // We’ll add all edges in random pair order
        let pairs = vec![(2,4), (4,1), (1,0), (0,3), (2,1), (1,3), (4,0), (4,3), (2,3), (0,1)];
        for (a, b) in pairs {
            // find the NodeIndex in map2
            let na = map2.iter().find(|(val, _)| *val == a).unwrap().1;
            let nb = map2.iter().find(|(val, _)| *val == b).unwrap().1;
            // If not already present, add edge
            if g2.find_edge(na, nb).is_none() {
                g2.add_edge(na, nb, ());
            }
        }

        // Now assert the two are isomorphic under canonization
        assert_isomorphic_canon(vec![g1, g2]);
    }

    #[test]
    fn test_bipartite_k33() {
        // K_{3,3} is the complete bipartite graph with partitions of size 3 each,
        // another symmetrical structure. We'll build two permutations and check isomorphism.

        let mut g1 = UnGraph::<char, ()>::new_undirected();
        // Partition A: a, b, c
        // Partition B: x, y, z
        let a = g1.add_node('a');
        let b = g1.add_node('b');
        let c = g1.add_node('c');
        let x = g1.add_node('x');
        let y = g1.add_node('y');
        let z = g1.add_node('z');
        // Connect all cross-partition edges
        for &left in &[a, b, c] {
            for &right in &[x, y, z] {
                g1.add_edge(left, right, ());
            }
        }

        // Another insertion order
        let mut g2 = UnGraph::<char, ()>::new_undirected();
        let z2 = g2.add_node('z');
        let b2 = g2.add_node('b');
        let a2 = g2.add_node('a');
        let x2 = g2.add_node('x');
        let c2 = g2.add_node('c');
        let y2 = g2.add_node('y');
        // Similarly connect cross-partitions but in random order:
        // The sets: {b2, a2, c2} and {z2, x2, y2}
        let edges = vec![(b2,x2), (b2,y2), (b2,z2),
                         (a2,x2), (a2,z2), (a2,y2),
                         (c2,x2), (c2,y2), (c2,z2)];
        for (m, n) in edges {
            g2.add_edge(m, n, ());
        }

        assert_isomorphic_canon(vec![g1, g2]);
    }

    #[test]
    fn test_disconnected_subgraphs() {
        // Two separate components, each a small cycle. The labeling algorithm must
        // handle disconnected graphs correctly. We'll build isomorphic variants
        // with different permutations of the components.

        // Graph 1: Two squares (4-cycles) not connected to each other.
        let mut g1 = UnGraph::<char, ()>::new_undirected();
        // First square: A-B-C-D
        let a1 = g1.add_node('A');
        let b1 = g1.add_node('B');
        let c1 = g1.add_node('C');
        let d1 = g1.add_node('D');
        g1.add_edge(a1, b1, ());
        g1.add_edge(b1, c1, ());
        g1.add_edge(c1, d1, ());
        g1.add_edge(d1, a1, ());
        // Second square: E-F-G-H
        let e1 = g1.add_node('E');
        let f1 = g1.add_node('F');
        let g11 = g1.add_node('G');
        let h1 = g1.add_node('H');
        g1.add_edge(e1, f1, ());
        g1.add_edge(f1, g11, ());
        g1.add_edge(g11, h1, ());
        g1.add_edge(h1, e1, ());

        // Graph 2: Same structure, different insertion order
        let mut g2 = UnGraph::<char, ()>::new_undirected();
        // We'll create the second square first
        let e2 = g2.add_node('E');
        let f2 = g2.add_node('F');
        let g2g = g2.add_node('G');
        let h2 = g2.add_node('H');
        g2.add_edge(e2, f2, ());
        g2.add_edge(f2, g2g, ());
        g2.add_edge(g2g, h2, ());
        g2.add_edge(h2, e2, ());
        // Then the first square
        let a2 = g2.add_node('A');
        let d2 = g2.add_node('D');
        let b2 = g2.add_node('B');
        let c2 = g2.add_node('C');
        g2.add_edge(a2, d2, ());
        g2.add_edge(d2, c2, ());
        g2.add_edge(c2, b2, ());
        g2.add_edge(b2, a2, ());

        assert_isomorphic_canon(vec![g1, g2]);
    }

    #[test]
    fn test_larger_cycle_of_size_8() {
        // A cycle of 8 nodes is highly symmetrical; ensuring that the canonical
        // labeling is stable can be tricky. We'll reorder the creation of nodes
        // and edges in multiple ways to see if we still get the same final canonical form.

        // Graph 1: straightforward 8-cycle (0 -> 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> back to 0)
        let mut g1 = UnGraph::<usize, ()>::new_undirected();
        let n: Vec<_> = (0..8).map(|i| g1.add_node(i)).collect();
        for i in 0..8 {
            g1.add_edge(n[i], n[(i + 1) % 8], ());
        }

        // Graph 2: same cycle but random insertion order
        let mut g2 = UnGraph::<usize, ()>::new_undirected();
        let indices = vec![4,1,0,7,3,6,2,5];
        let nodes2: Vec<_> = indices.iter().map(|&i| g2.add_node(i)).collect();
        // connect them in a cycle order based on the actual integer label
        // (sort-of “random” but consistent)
        // We'll connect i->(i+1) mod 8 in terms of *labels*, not node order
        // so we need a function to get NodeIndex by label:
        fn find_node_by_label(g: &UnGraph<usize, ()>, label: usize) -> NodeIndex {
            g.node_indices()
                .find(|idx| *g.node_weight(*idx).unwrap() == label)
                .unwrap()
        }

        for i in 0..8 {
            let a_label = i;
            let b_label = (i + 1) % 8;
            let na = find_node_by_label(&g2, a_label);
            let nb = find_node_by_label(&g2, b_label);
            if g2.find_edge(na, nb).is_none() {
                g2.add_edge(na, nb, ());
            }
        }

        assert_isomorphic_canon(vec![g1, g2]);
    }

    #[test]
    fn test_multiple_isomorphic_subgraphs() {
        // Build a graph that has repeated isomorphic “subgraphs” within the same graph:
        // i.e., 3 disconnected components, all triangles. Then reorder how we add them,
        // ensuring canonical labeling does not break.

        let mut g1 = UnGraph::<char, ()>::new_undirected();
        // tri1: A-B, B-C, A-C
        let a1 = g1.add_node('A');
        let b1 = g1.add_node('B');
        let c1 = g1.add_node('C');
        g1.add_edge(a1, b1, ());
        g1.add_edge(b1, c1, ());
        g1.add_edge(c1, a1, ());
        // tri2: D-E-F
        let d1 = g1.add_node('D');
        let e1 = g1.add_node('E');
        let f1 = g1.add_node('F');
        g1.add_edge(d1, e1, ());
        g1.add_edge(e1, f1, ());
        g1.add_edge(f1, d1, ());
        // tri3: X-Y-Z
        let x1 = g1.add_node('X');
        let y1 = g1.add_node('Y');
        let z1 = g1.add_node('Z');
        g1.add_edge(x1, y1, ());
        g1.add_edge(y1, z1, ());
        g1.add_edge(z1, x1, ());

        // Another variant
        let mut g2 = UnGraph::<char, ()>::new_undirected();
        // tri2 first
        let e2 = g2.add_node('E');
        let d2 = g2.add_node('D');
        let f2 = g2.add_node('F');
        g2.add_edge(e2, f2, ());
        g2.add_edge(f2, d2, ());
        g2.add_edge(d2, e2, ());
        // tri3 next
        let z2 = g2.add_node('Z');
        let x2 = g2.add_node('X');
        let y2 = g2.add_node('Y');
        g2.add_edge(z2, x2, ());
        g2.add_edge(x2, y2, ());
        g2.add_edge(y2, z2, ());
        // tri1 last
        let b2 = g2.add_node('B');
        let c2 = g2.add_node('C');
        let a2 = g2.add_node('A');
        g2.add_edge(b2, a2, ());
        g2.add_edge(a2, c2, ());
        g2.add_edge(c2, b2, ());

        assert_isomorphic_canon(vec![g1, g2]);
    }

    //
    // If desired, you can add tests with more nodes & edges or random graphs.
    // However, random tests should fix a seed to remain deterministic for unit tests.
    //
}
 */
