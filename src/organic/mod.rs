use std::{collections::VecDeque, fmt::{Debug, Formatter, Result as FmtResult}, ops::RangeBounds};

use petgraph::graph::Node;
use anyhow::{Result, Context, anyhow, bail};
use super::*;
use tracing::*;

fn remove_subdivision<N, E>(graph: &mut UnGraph<N, E>, node: NodeIndex) where E: PartialEq + Clone {
    // Get neighbors
    let neighbors: Vec<NodeIndex> = graph.neighbors(node).collect();

    // Ensure the node is a subdivision (has exactly two neighbors)
    if neighbors.len() == 2 {
        // Get the weight of the edge between the neighbors
        let edge1 = graph.find_edge(node, neighbors[0]).unwrap();
        let edge2 = graph.find_edge(node, neighbors[1]).unwrap();
        let are_equal = graph[edge1] == graph[edge2];
        let weight = graph[edge1].clone();
        // Connect the two neighbors
        graph.add_edge(neighbors[0], neighbors[1], weight);

        graph.remove_edge(edge1).unwrap();
        graph.remove_edge(edge2).unwrap();

        if !are_equal {
            panic!("Removing subdivision between inconsistent edge weights");
        }

        // Remove the subdivision node
        graph.remove_node(node);
    }
}

#[derive(Clone)]
pub struct Substituent(MoleculeGraph);

impl Substituent {
    pub fn is_same_as(&self, other: &Substituent) -> bool {
        petgraph::algo::is_isomorphic(&self.0, &other.0)
    }

    pub fn to_molecule_graph(&self) -> MoleculeGraph {
        // First, remove any R groups that are connected
        let mut clone = self.clone();
        clone.remove_connected_r_groups();
        clone.0
    }

    pub fn remove_all_r_groups(mut self) -> Self {
        self.0.retain_nodes(|g, node| !g[node].is_r_group());
        self
    }

    fn to_raw_graph(&self) -> &MoleculeGraph {
        &self.0
    }

    pub fn visualize_raw(&self, filename: &str) -> Result<()> {
        let mut graph = self.to_raw_graph().clone();
        graph.hydrogenate();
        visualize_graph(&graph, "tmp.dot", Some(filename)).map_err(|e| anyhow!("Failed to visualize: {}", e))?;
        Ok(())
    }

    pub fn visualize(&self, filename: &str) -> Result<()> {
        let mut graph = self.to_molecule_graph();
        graph.hydrogenate();
        visualize_graph(&graph, "tmp.dot", Some(filename)).map_err(|e| anyhow!("Failed to visualize: {}", e))?;
        Ok(())
    }
    /// Parse a SMILES string into a substituent.
    ///
    /// This minimal parser supports:
    /// - Atoms (e.g. C, O, N, Cl, Br)
    /// - The literal `R` for R groups (parsed as an R group node)
    /// - Optional explicit bonds: '-' (single), '=' (double), '#' (triple)
    /// - Branching via parentheses
    /// - Ring closures using digits
    pub fn parse_smiles(smiles: &str) -> Result<Self> {
        use std::iter::Peekable;
        use std::str::Chars;

        // Create a new molecule graph.
        let mut graph = MoleculeGraph::default();

        // Mapping from ring closure digit to node index.
        let mut ring_closures: HashMap<char, NodeIndex> = HashMap::new();
        // A stack to hold nodes for branch starting points.
        let mut branch_stack: Vec<NodeIndex> = Vec::new();
        // If a bond is explicitly specified, we store it here.
        let mut pending_bond: Option<Bond> = None;
        // The most-recently parsed atom.
        let mut current_node: Option<NodeIndex> = None;

        // Create a peekable iterator over characters.
        let mut chars: Peekable<Chars> = smiles.chars().peekable();

        while let Some(ch) = chars.next() {
            match ch {
                // Start a branch: push the current atom onto the stack.
                '(' => {
                    if let Some(node) = current_node {
                        branch_stack.push(node);
                    } else {
                        return Err(anyhow!("Encountered '(' with no current atom"));
                    }
                }
                // End a branch: pop a saved node.
                ')' => {
                    current_node = Some(branch_stack.pop()
                        .ok_or_else(|| anyhow!("Unmatched ')' in SMILES"))?);
                }
                // Explicit bond symbols.
                '-' | '=' | '#' | ':' => {
                    pending_bond = Some(match ch {
                        '-' => Bond::Single,
                        '=' => Bond::Double,
                        '#' => Bond::Triple,
                        ':' => Bond::Aromatic,
                        _ => unreachable!(),
                    });
                }
                // A digit indicates a ring closure.
                d if d.is_ascii_digit() => {
                    if let Some(stored) = ring_closures.remove(&d) {
                        // We have seen this digit before—close the ring.
                        if let Some(curr) = current_node {
                            let bond = pending_bond.take().unwrap_or(Bond::Single);
                            graph.add_edge(stored, curr, bond);
                        } else {
                            return Err(anyhow!("Ring closure digit encountered with no current atom"));
                        }
                    } else {
                        // First occurrence: store the current atom for later connection.
                        if let Some(curr) = current_node {
                            ring_closures.insert(d, curr);
                        } else {
                            return Err(anyhow!("Ring closure digit encountered with no current atom"));
                        }
                    }
                }
                // An alphabetic character starts an atom specification.
                ch if ch.is_alphabetic() => {
                    // Special case: the literal "R" represents an R group.
                    if ch == 'R' {
                        // Count the following number of apostrophes to determine the R group number.
                        let mut count = 0;
                        while let Some(&'\'') = chars.peek() {
                            count += 1;
                            chars.next();
                        }
                        debug!("R group count: {}", count);
                        let r_group = Element::r_group(count);
                        let new_node = graph.add_node(r_group);
                        if let Some(prev) = current_node {
                            let bond = pending_bond.take().unwrap_or(Bond::Single);
                            graph.add_edge(prev, new_node, bond);
                        }
                        current_node = Some(new_node);
                    } else {
                        // Parse an element symbol. If the next character is lowercase, include it.
                        let mut symbol = ch.to_string();
                        if ch.is_ascii_uppercase() {
                            if let Some(&next) = chars.peek() {
                                if next.is_lowercase() {
                                    symbol.push(chars.next().unwrap());
                                }
                            }
                        }

                        let element = match symbol.as_str() {
                            "C"  => Element::C,
                            "O"  => Element::O,
                            "N"  => Element::N,
                            "Cl" => Element::Cl,
                            "Br" => Element::Br,
                            "c"  => Element::C.as_aromatic(),
                            "o"  => Element::O.as_aromatic(),
                            "n"  => Element::N.as_aromatic(),
                            // Extend here for additional elements as needed.
                            _ => return Err(anyhow!("Unknown element symbol: {}", symbol)),
                        };
                        let new_node = graph.add_node(element);
                        if let Some(prev) = current_node {
                            let bond = pending_bond.take().unwrap_or(Bond::Single);
                            graph.add_edge(prev, new_node, bond);
                        }
                        current_node = Some(new_node);
                    }
                }
                // Anything else is unexpected.
                _ => {
                    return Err(anyhow!("Unexpected character in SMILES: {}", ch));
                }
            }
        }

        if !branch_stack.is_empty() {
            return Err(anyhow!("Unmatched '(' found in SMILES"));
        }
        if !ring_closures.is_empty() {
            return Err(anyhow!("Unclosed ring closures in SMILES: {:?}", ring_closures.keys()));
        }

        Ok(Substituent(graph))
    }

    /// Remove any R groups that are connected -- these are just placeholders for
    /// connecting substituents together that we eliminate after connecting them.
    fn remove_connected_r_groups(&mut self) {
        // Find all the R groups
        while self.has_connected_r_groups() {
            for r_group in self.find_r_groups() {
                // Check to confirm that the R group is connected to more than one node
                if self.0.neighbors(r_group).count() <= 1 {
                    continue;
                }
                remove_subdivision(&mut self.0, r_group);
                break
            }
        }
    }

    /// Does this have any R groups that we need to remove?
    fn has_connected_r_groups(&self) -> bool {
        self.find_r_groups().iter().any(|&r_group| self.0.neighbors(r_group).count() > 1)
    }

    /// Find all the R groups in the substituent
    fn find_r_groups(&self) -> Vec<NodeIndex> {
        // Find all the node indices that are R groups
        let mut node_indices: Vec<NodeIndex> = self.0.node_indices().filter(|&node| self.0[node].is_r_group()).collect();

        // Sort them by the R group number
        node_indices.sort_by_key(|&node| match self.0[node].kind {
            ElementType::RGroup(n) => n,
            _ => unreachable!(),
        });

        node_indices
    }

    /// Find the lowest R group in the substituent (the next one to connect)
    fn find_lowest_r_group(&self) -> Option<NodeIndex> {
        self.find_r_groups().first().copied()
    }

    fn duplicate_r_group(&mut self, r_group: NodeIndex) -> NodeIndex {
        // Find the element the R group is connected to
        let element = self.0.neighbors(r_group).find(|&node| !self.0[node].is_r_group()).unwrap();
        let new_r_group = self.0.add_node(Element::r_group(self.count_r_groups()));
        self.0.add_edge(element, new_r_group, Bond::Single);
        new_r_group
    }

    /// Does this have any R groups?
    fn has_r_groups(&self) -> bool {
        !self.find_r_groups().is_empty()
    }

    fn count_r_groups(&self) -> usize {
        self.find_r_groups().len()
    }

    fn connect_atom(&mut self, atom: Element, at: NodeIndex) {
        let new_atom = self.0.add_node(atom);
        self.0.add_edge(at, new_atom, Bond::Single);
    }
    
    /// Connect another substituent to this one.
    /// 
    /// This will connext the next lowest R group in this substituent to the next lowest R group in the other substituent.
    fn connect(&mut self, other: &Substituent) -> Result<()> {
        if self.count_r_groups() > 1 && other.count_r_groups() > 1 {
            warn!("Connecting two substituents with multiple R groups, future connections may be in an unintended order!");
        }

        // Find the R group in the other substituent
        let my_r_group = self.find_lowest_r_group().context("No R groups found")?;
        let other_r_group = other.find_lowest_r_group().context("No R groups found")?;
        super::connect(&mut self.0, &other.0, my_r_group, other_r_group, Bond::Single);
        self.remove_connected_r_groups();
        Ok(())
    }

    /// Connect another substituent to this one at a specific atom.
    /// 
    /// This uses the R group in the other substituent and connects it to the specified atom in this substituent.
    fn connect_at(&mut self, other: &Substituent, my_atom: NodeIndex) -> Result<()> {
        let other_r_group = other.find_lowest_r_group().context("No R groups found")?;
        super::connect(&mut self.0, &other.0, my_atom, other_r_group, Bond::Single);
        self.remove_connected_r_groups();
        Ok(())
    }

    pub fn then_connect(mut self, other: Substituent) -> Result<Self> {
        assert!(self.has_r_groups(), "No R groups found in the first substituent");
        assert!(other.has_r_groups(), "No R groups found in the second substituent");

        self.connect(&other)?;
        Ok(self)
    }

    /// Find the longest carbon chain in the substituent.
    pub fn find_longest_carbon_chain(&self) -> Vec<NodeIndex> {
        // Perform a DFS to find the longest carbon chain
        let mut dfs = petgraph::visit::Dfs::new(&self.0, self.0.node_indices().next().unwrap());
        let mut longest_chain = Vec::new();
        let mut current_chain = Vec::new();
        let mut last_node = None;
        while let Some(node) = dfs.next(&self.0) {
            if let Some(last_node) = last_node {
                // Confirm that the current node is connected to the last node, else reset the chain
                if !self.0.neighbors(last_node).any(|n| n == node) {
                    current_chain.clear();
                }
            }

            if self.0[node].kind == ElementType::C {
                current_chain.push(node);
            } else {
                if current_chain.len() > longest_chain.len() {
                    longest_chain = current_chain.clone();
                }
                current_chain.clear();
            }

            last_node = Some(node);
        }

        if current_chain.len() > longest_chain.len() {
            longest_chain = current_chain;
        }

        longest_chain
    }
    
    /// Is this a linear chain?
    pub fn is_chain(&self) -> bool {
        if self.is_cyclic() {
            return false;
        }

        if petgraph::algo::connected_components(&self.0) > 1 {
            return false;
        }

        // Confirm this is a linear chain -- each node has at most two neighbors, except for the ends
        let mut ends = 0;
        for node in self.0.node_indices() {
            let neighbors: Vec<NodeIndex> = self.0.neighbors(node).collect();
            match neighbors.len() {
                1 => ends += 1,
                2 => continue,
                _ => return false,
            }
        }

        ends == 2
    }

    /// Is this a ring?
    pub fn is_ring(&self) -> bool {
        if !self.is_cyclic() || petgraph::algo::connected_components(&self.0) != 1 {
            return false;
        }

        // Confirm this is a ring -- each node has exactly two neighbors
        for node in self.0.node_indices() {
            let neighbors: Vec<NodeIndex> = self.0.neighbors(node).collect();
            if neighbors.len() != 2 {
                return false;
            }
        }

        true
    }

    /// Is this cyclic?
    pub fn is_cyclic(&self) -> bool {
        petgraph::algo::is_cyclic_undirected(&self.0)
    }


    /*
    /// Create a new substituent with N carbons in a row all connected by single bonds.
    pub fn alkane(n: u32) -> Self {
        let mut graph = MoleculeGraph::default();

        let carbon = graph.add_node(Element::C);
        let mut prev = carbon;

        for _ in 1..n {
            let next = graph.add_node(Element::C);
            graph.add_edge(prev, next, Bond::Single);
            prev = next;
        }

        Substituent(graph)
    }

    /// Create a new substituent with N carbons in a ring all connected by aromatic bonds.
    pub fn cyclic_alkane(n: u32) -> Self {
        let mut graph = MoleculeGraph::default();

        let carbon = graph.add_node(Element::C);
        let mut prev = carbon;

        let bond = if n == 6 {
            Bond::Aromatic
        } else {
            Bond::Single
        };

        for _ in 1..n {
            let next = graph.add_node(Element::C);
            graph.add_edge(prev, next, bond);
            prev = next;
        }

        graph.add_edge(prev, carbon, bond);

        Substituent(graph)
    }
     */

    pub fn alkane(n: usize) -> Self {
        Self::alkyl_group(n, vec![])
    }

    pub fn alkyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        Self::carbon_chain_group_helper(n, vec![Bond::Single; n - 1], r_group_positions, false)
    }

    pub fn alkenyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        Self::carbon_chain_group_helper(n, vec![Bond::Double; n - 1], r_group_positions, false)
    }

    pub fn alkynyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        Self::carbon_chain_group_helper(n, vec![Bond::Triple; n - 1], r_group_positions, false)
    }

    fn carbon_chain_group_helper(n: usize, bonds: impl IntoIterator<Item=Bond>, r_group_positions: impl IntoIterator<Item = usize>, cyclic: bool) -> Self {
        if cyclic && n < 3 {
            panic!("Cyclic alkyl groups must have at least 3 carbons");
        }

        let mut graph = MoleculeGraph::default();
        let r_group_positions: Vec<usize> = r_group_positions.into_iter().collect();
        let bonds: Vec<Bond> = bonds.into_iter().collect();
        let mut prev = graph.add_node(Element::C);
        let mut count = 0;

        let mut nodes = vec![prev];

        for i in 1..n {
            let next = graph.add_node(Element::C);
            graph.add_edge(prev, next, bonds[i-1]);
            nodes.push(next);

            prev = next;
        }

        if cyclic {
            graph.add_edge(prev, nodes[0], bonds[n-1]);
        }

        for pos in r_group_positions {
            if pos > n {
                panic!("R group position {} is greater than the number of carbons in the chain", pos);
            }
            let node = nodes[pos - 1];
            let r_group = graph.add_node(Element::r_group(count));
            count += 1;
            graph.add_edge(node, r_group, Bond::Single);
        }

        Substituent(graph)
    }

    pub fn phenyl_group() -> Self {
        Self::cyclic_alkyl_group(6, vec![1])
    }

    pub fn carbonyl_group() -> Self {
        let mut graph = MoleculeGraph::default();
        let carbon = graph.add_node(Element::C);
        let oxygen = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen, Bond::Double);

        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(carbon, r1, Bond::Single);
        let r2 = graph.add_node(Element::r_group(1));
        graph.add_edge(carbon, r2, Bond::Single);

        Substituent(graph)
    }

    pub fn carboxyl_group() -> Self {
        let mut graph = MoleculeGraph::default();
        let carbon = graph.add_node(Element::C);
        let oxygen1 = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen1, Bond::Double);
        let oxygen2 = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen2, Bond::Single);

        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(carbon, r1, Bond::Single);

        Substituent(graph)
    }

    pub fn ketone() -> Self {
        let mut graph = MoleculeGraph::default();
        let carbon = graph.add_node(Element::C);
        let oxygen = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen, Bond::Double);

        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(carbon, r1, Bond::Single);

        let r2 = graph.add_node(Element::r_group(1));
        graph.add_edge(carbon, r2, Bond::Single);

        Substituent(graph)
    }

    pub fn aldehyde() -> Self {
        let mut graph = MoleculeGraph::default();
        let carbon = graph.add_node(Element::C);
        let oxygen = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen, Bond::Double);

        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(carbon, r1, Bond::Single);

        Substituent(graph)
    }

    pub fn ester() -> Self {
        let mut graph = MoleculeGraph::default();
        let carbon = graph.add_node(Element::C);
        let oxygen1 = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen1, Bond::Double);
        let oxygen2 = graph.add_node(Element::O);
        graph.add_edge(carbon, oxygen2, Bond::Single);

        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(carbon, r1, Bond::Single);
        let r2 = graph.add_node(Element::r_group(1));
        graph.add_edge(oxygen2, r2, Bond::Single);

        Substituent(graph)
    }

    pub fn ether() -> Self {
        let mut graph = MoleculeGraph::default();
        let oxygen = graph.add_node(Element::O);
        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(oxygen, r1, Bond::Single);
        let r2 = graph.add_node(Element::r_group(1));
        graph.add_edge(oxygen, r2, Bond::Single);

        Substituent(graph)
    }

    pub fn cyclic_alkyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        Self::carbon_chain_group_helper(n, vec![Bond::Single; n], r_group_positions, true)
    }

    pub fn cyclic_alkane(n: usize) -> Self {
        Self::carbon_chain_group_helper(n, vec![Bond::Single; n], vec![], true)
    }
    
    pub fn methane() -> Self {
        Self::alkane(1)
    }

    pub fn ethane() -> Self {
        Self::alkane(2)
    }

    pub fn propane() -> Self {
        Self::alkane(3)
    }

    pub fn butane() -> Self {
        Self::alkane(4)
    }

    pub fn pentane() -> Self {
        Self::alkane(5)
    }

    pub fn hexane() -> Self {
        Self::alkane(6)
    }

    pub fn methyl() -> Self {
        Self::alkyl_group(1, vec![1])
    }

    pub fn ethyl() -> Self {
        Self::alkyl_group(2, vec![1])
    }

    pub fn propyl() -> Self {
        Self::alkyl_group(3, vec![1])
    }

    pub fn butyl() -> Self {
        Self::alkyl_group(4, vec![1])
    }

    pub fn pentyl() -> Self {
        Self::alkyl_group(5, vec![1])
    }

    pub fn hexyl() -> Self {
        Self::alkyl_group(6, vec![1])
    }

    pub fn isopropyl() -> Self {
        Self::alkyl_group(3, vec![2])
    }

    pub fn isobutyl() -> Self {
        Substituent::alkyl_group(3, vec![2])
            .then_connect(Substituent::alkyl_group(1, vec![1, 1])).unwrap()
    }

    pub fn secbutyl() -> Self {
        Substituent::alkyl_group(3, vec![1, 2])
            .then_connect(Substituent::methyl()).unwrap()
    }

    pub fn tertbutyl() -> Self {
        Substituent::alkyl_group(1, vec![1, 1, 1, 1])
            .then_connect(Substituent::methyl()).unwrap()
            .then_connect(Substituent::methyl()).unwrap()
            .then_connect(Substituent::methyl()).unwrap()
    }

    pub fn halide(elem: impl Into<ElementType>) -> Self {
        let elem = elem.into();
        let mut graph = MoleculeGraph::default();
        let halide = graph.add_node(elem.into());
        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(r1, halide, Bond::Single);
        Substituent(graph)
    }

    pub fn hydrogen() -> Self {
        let mut graph = MoleculeGraph::default();
        let halide = graph.add_node(Element::H);
        let r1 = graph.add_node(Element::r_group(0));
        graph.add_edge(r1, halide, Bond::Single);
        Substituent(graph)
    }
}

impl From<OrganicMolecule> for Substituent {
    fn from(mol: OrganicMolecule) -> Self {
        mol.as_substituent()
    }
}

impl From<Element> for Substituent {
    fn from(elem: Element) -> Self {
        let mut graph = MoleculeGraph::default();
        let node = graph.add_node(elem);
        let r = graph.add_node(Element::r_group(0));
        graph.add_edge(node, r, Bond::Single);
        Substituent(graph)
    }
}

impl Debug for Substituent {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", molecule_to_smiles(&self.0))
    }
}


/// --- Helper Functions ---
///
/// The following functions are sketches for operations on MoleculeGraph. You would need
/// to implement them (or adapt petgraph functions) for your actual data structures.

/// Splits the given graph by removing the edge with the specified index and
/// returns the two resulting connected subgraphs.
fn split_graph_at_edge(graph: &MoleculeGraph, edge: EdgeIndex) -> (MoleculeGraph, MoleculeGraph) {
    // Clone the graph so we can modify it.
    let mut g = graph.clone();
    // Get the nodes connected by the edge.
    let (source, target) = graph.edge_endpoints(edge).unwrap();

    let r1 = g.add_node(Element::r_group(0));
    let r2 = g.add_node(Element::r_group(0));
    g.add_edge(source, r1, Bond::Single);
    g.add_edge(target, r2, Bond::Single);
    g.remove_edge(edge).expect("Edge removal failed");

    // Use petgraph’s connected components (or any other suitable algorithm) to partition the graph.
    // For simplicity, assume there are exactly two components.
    let comps = petgraph::algo::connected_components(&g);

    // (Pseudo‑code) Extract subgraphs corresponding to each component.
    let subgraph1 = extract_subgraph(&g, /*component_id=*/ 0);
    let subgraph2 = extract_subgraph(&g, /*component_id=*/ 1);

    // Add R groups to the subgraphs
    // let r1 = subgraph1.add_node(Element::r_group(0));
    // let r2 = subgraph2.add_node(Element::r_group(1));

    (subgraph1, subgraph2)
}

/// Removes the specified node from the graph and returns a new graph.
fn remove_node_from_graph(graph: &MoleculeGraph, node: NodeIndex) -> MoleculeGraph {
    let mut g = graph.clone();
    g.remove_node(node);
    g
}

/// Returns a vector of components (each represented by a set of NodeIndex) for the graph.
fn connected_components(graph: &MoleculeGraph) -> Vec<Vec<NodeIndex>> {
    // We'll perform a DFS starting at each unvisited node.
    let mut components: Vec<Vec<NodeIndex>> = Vec::new();
    let mut visited: HashSet<NodeIndex> = HashSet::new();

    for node in graph.node_indices() {
        if visited.contains(&node) {
            continue;
        }
        // Start a new component.
        let mut stack = vec![node];
        let mut comp = Vec::new();
        while let Some(current) = stack.pop() {
            if visited.insert(current) {
                comp.push(current);
                for neighbor in graph.neighbors(current) {
                    if !visited.contains(&neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
        components.push(comp);
    }
    components
}

/// Given a graph and a set of node indices representing one component, returns a subgraph.
///
/// The new graph will contain copies of only the nodes in `component` and any edge
/// that connected two nodes in `component`. Note that the node indices in the subgraph
/// will not match those in the original graph.
fn subgraph_from_component(graph: &MoleculeGraph, component: &Vec<NodeIndex>) -> MoleculeGraph {
    let mut subgraph = MoleculeGraph::default();
    // Mapping from old NodeIndex to new NodeIndex in the subgraph.
    let mut node_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();

    // Insert all nodes from the component.
    for &node in component.iter() {
        // Assume the node weight (e.g. an Element) implements Clone.
        let new_node = subgraph.add_node(graph[node].clone());
        node_map.insert(node, new_node);
    }

    // For each edge in the original graph, if both endpoints are in the component,
    // add the corresponding edge to the subgraph.
    for edge in graph.edge_references() {
        let source = edge.source();
        let target = edge.target();
        if component.contains(&source) && component.contains(&target) {
            let new_source = node_map[&source];
            let new_target = node_map[&target];
            // Assume Bond (the edge weight) implements Clone.
            subgraph.add_edge(new_source, new_target, edge.weight().clone());
        }
    }
    subgraph
}

/// Extract a subgraph from the graph corresponding to the given component ID.
/// This function computes the connected components and returns the subgraph
/// for the component with index `_component_id`.
fn extract_subgraph(graph: &MoleculeGraph, _component_id: usize) -> MoleculeGraph {
    let comps = connected_components(graph);
    if _component_id >= comps.len() {
        panic!(
            "Component id {} is out of range; found only {} components",
            _component_id,
            comps.len()
        );
    }
    subgraph_from_component(graph, &comps[_component_id])
}

/// Finds a subgraph isomorphism from `pattern` (a MoleculeGraph) into `graph`.
/// Returns a mapping from each node in `pattern` to its corresponding node in `graph`.
fn find_subgraph_isomorphism(
    graph: &MoleculeGraph,
    pattern: &MoleculeGraph,
) -> Option<HashMap<NodeIndex, NodeIndex>> {
    // Collect pattern nodes in a vector.
    let pattern_nodes: Vec<NodeIndex> = pattern.node_indices().collect();
    let mut mapping = HashMap::new();
    let mut used = HashSet::new();

    fn backtrack(
        graph: &MoleculeGraph,
        pattern: &MoleculeGraph,
        pattern_nodes: &[NodeIndex],
        mapping: &mut HashMap<NodeIndex, NodeIndex>,
        used: &mut HashSet<NodeIndex>,
        index: usize,
    ) -> bool {
        if index == pattern_nodes.len() {
            return true;
        }
        let p_node = pattern_nodes[index];
        for g_node in graph.node_indices() {
            if used.contains(&g_node) {
                continue;
            }
            // // Check node compatibility (here we simply require the same Element type).
            // if graph[g_node].kind != pattern[p_node].kind {
            //     continue;
            // }

            // If the pattern node is not an R group, enforce that the element types match.
            if !pattern[p_node].is_r_group() && graph[g_node].kind != pattern[p_node].kind {
                continue;
            }
            // If the pattern node is an R group, enforce that the matched atom is not an R group.
            if pattern[p_node].is_r_group() && graph[g_node].is_r_group() {
                continue;
            }
            // For every neighbor of p_node already mapped, verify that the corresponding edge exists.
            let mut valid = true;
            for edge in pattern.edges(p_node) {
                let p_nbr = edge.target();
                if let Some(&mapped) = mapping.get(&p_nbr) {
                    let found = graph.edges(g_node)
                        .any(|e| e.target() == mapped && *e.weight() == *edge.weight());
                    if !found {
                        valid = false;
                        break;
                    }
                }
            }
            if !valid {
                continue;
            }
            mapping.insert(p_node, g_node);
            used.insert(g_node);
            if backtrack(graph, pattern, pattern_nodes, mapping, used, index + 1) {
                return true;
            }
            mapping.remove(&p_node);
            used.remove(&g_node);
        }
        false
    }

    if backtrack(graph, pattern, &pattern_nodes, &mut mapping, &mut used, 0) {
        // Remove mapping entries for R group nodes (wildcards).
        mapping.retain(|&p, _| !pattern[p].is_r_group());
        Some(mapping)
    } else {
        None
    }
}

/// Disconnects from `graph` any subgraph isomorphic to the given `pattern` (passed as a Substituent).
/// It “cuts off” all edges between the matched nodes and the rest of the graph,
/// then returns the connected components (each wrapped as a Substituent).
fn disconnect_subgraph_by_pattern(
    graph: &MoleculeGraph,
    pattern_sub: &Substituent,
) -> Result<(Vec<Substituent>, Substituent)> {
    // Remove all R groups that are connected to the pattern
    // let pattern_sub = pattern_sub.clone().remove_all_r_groups();
    // Find the matching mapping.
    let pattern = pattern_sub.to_raw_graph();

    let mapping = find_subgraph_isomorphism(graph, pattern)
        .ok_or_else(|| anyhow!("No matching subgraph found for the given pattern"))?;
    let matched_set: HashSet<NodeIndex> = mapping.values().copied().collect();
    let inverted_mapping: HashMap<NodeIndex, NodeIndex> = mapping.into_iter().map(|(k, v)| (v, k)).collect();
    
    // Clone the graph so we can modify it.
    let mut g = graph.clone();
    // Remove all edges that connect a matched node to a non‐matched node.
    let edges_to_remove: Vec<EdgeIndex> = g.edge_references()
        .filter(|e| {
            let src = e.source();
            let tgt = e.target();
            (matched_set.contains(&src) && !matched_set.contains(&tgt))
                || (matched_set.contains(&tgt) && !matched_set.contains(&src))
        })
        .map(|e| e.id())
        .collect();

    let cut_nodes: HashSet<NodeIndex> = edges_to_remove.iter()
        .flat_map(|&edge| {
            let (src, tgt) = g.edge_endpoints(edge).unwrap();
            vec![src, tgt]
        })
        .filter(|&node| !matched_set.contains(&node))
        .collect();

    g.retain_edges(|_, edge| !edges_to_remove.contains(&edge));


    // (We leave the matched nodes in the graph here.)
    // Return the connected components.
    let mut comps = connected_components(&g);

    // Sort the connected components by which node they were connected to
    comps.sort_by_key(|comp| {
        if comp.is_empty() {
            return -1;
        }

        if comp.iter().any(|n| matched_set.contains(n)) {
            return 0;
        }

        // Find the node in the matched set that the component is connected to
        if let Some(connected_node) = graph.node_indices().find(|&node| {
            // If this node is in the matched set, and its connected to the component
            matched_set.contains(&node) && graph.neighbors(node).any(|n| comp.contains(&n))
        }) {
            // Find the index of the connected node in the matched set
            inverted_mapping[&connected_node].index() as i32
        } else {
            panic!("Component not connected to any matched nodes");
        }
    });

    let mut component_r_group_count = HashMap::new();
    for node in cut_nodes {
        // Get the connected component
        let comp = comps.iter_mut().find(|comp| comp.contains(&node)).unwrap();
        // Count the number of R groups in the component
        let r_group_count = comp.iter().filter(|&n| g[*n].is_r_group()).count();
        // Store the count for the component
        component_r_group_count.insert(node, r_group_count);

        // Add an R group to the cut node
        let r_group = g.add_node(Element::r_group(r_group_count));
        g.add_edge(node, r_group, Bond::Single);

        // Add the R group to the component
        comp.push(r_group);
    }

    // let comps = connected_components(&g);
    let pattern_sub = comps
        .iter()
        .find(|comp| comp.iter().any(|n| matched_set.contains(n)))
        .map(|comp| Substituent(subgraph_from_component(&g, &comp)))
        .unwrap();
    
    let subs: Vec<Substituent> = comps.into_iter()
        .filter(|comp| {
            // Confirm that the component didn't contain any of the matched nodes
            comp.iter().all(|node| !matched_set.contains(node))
        })
        .map(|comp| Substituent(subgraph_from_component(&g, &comp)))
        .collect();

    Ok((subs, pattern_sub))
}


#[derive(Clone)]
pub enum OrganicMolecule {
    Substituent(Substituent),
    Compound(UnGraph<Self, usize>)
}

impl OrganicMolecule {
    pub fn is_same_as(&self, other: &Self) -> bool {
        self.as_substituent().is_same_as(&other.as_substituent())
    }

    /// Parse a SMILES string into an OrganicMolecule by first building its graph and then
    /// structurally matching known functional group patterns.
    pub fn parse_smiles(smiles: &str) -> Result<Self> {
        // First, parse the SMILES into a substituent (and hence a MoleculeGraph)
        let substituent = Substituent::parse_smiles(smiles)?;
        let raw_graph = substituent.to_raw_graph().clone();

        // Try to match known functional groups structurally.
        info!("Matching functional groups in the molecule");
        info!("Matching ester...");
        if let Some((acyl_graph, alkyl_graph)) = Self::match_ester(&raw_graph) {
            let acyl = OrganicMolecule::from(Substituent(acyl_graph));
            let alkyl = OrganicMolecule::from(Substituent(alkyl_graph));
            info!("Ester detected: acyl = {:?}, alkyl = {:?}", acyl, alkyl);
            return Ok(OrganicMolecule::ester(acyl, alkyl));
        }
        info!("Matching ether...");
        if let Some((alkyl1_graph, alkyl2_graph)) = Self::match_ether(&raw_graph) {
            let alkyl1 = OrganicMolecule::from(Substituent(alkyl1_graph));
            let alkyl2 = OrganicMolecule::from(Substituent(alkyl2_graph));
            info!("Ether detected: alkyl1 = {:?}, alkyl2 = {:?}", alkyl1, alkyl2);
            return Ok(OrganicMolecule::ether(alkyl1, alkyl2));
        }
        info!("Matching carboxylic acid...");
        if let Some(target_graph) = Self::match_carboxylic_acid(&raw_graph) {
            info!("Carboxylic acid detected: target = {:?}", target_graph);
            let target = OrganicMolecule::from(Substituent(target_graph));
            return Ok(OrganicMolecule::carboxylic_acid(target));
        }
        info!("No special functional groups detected");
        // If no special group was detected, fall back to a plain substituent.
        Ok(OrganicMolecule::from(substituent))
    }

    pub fn bind(&self, sites: impl IntoIterator<Item=Self>) -> Result<Self> {
        let sites = sites.into_iter().collect::<Vec<_>>();
        // Confirm we have enough sites to bind
        if self.count_r_groups() < sites.len() {
            bail!("Not enough R groups in the target molecule to bind all sites");
        }

        // Create a new organic molecule
        let mut graph = UnGraph::new_undirected();
        let root = graph.add_node(self.clone());

        // Bind the sites to the molecule
        for site in sites {
            let node = graph.add_node(site);
            graph.add_edge(root, node, 0);
        }

        Ok(OrganicMolecule::Compound(graph))
    }

    /// Match an ester functional group in the given graph.
    ///
    /// This method looks for a carbon atom that is bonded as follows:
    /// - Exactly one double-bonded oxygen neighbor (carbonyl oxygen)
    /// - Exactly one single-bonded oxygen neighbor (alkoxy oxygen)
    /// - Exactly one other neighbor (the acyl placeholder)
    ///
    /// In addition, the alkoxy oxygen must be bonded to exactly one other node,
    /// which becomes the alkyl placeholder.
    ///
    /// When found, the entire ester center (the carbon and both oxygens) is removed
    /// so that neither any part of the C(=O)-O structure remains in the fragments.
    /// The connected component containing the acyl placeholder becomes the acyl fragment,
    /// and that containing the alkyl placeholder becomes the alkyl fragment.
    fn match_ester(graph: &MoleculeGraph) -> Option<(MoleculeGraph, MoleculeGraph)> {
        let (subs, _) = disconnect_subgraph_by_pattern(graph, &Substituent::ester()).ok()?;
        if subs.len() == 2 {
            let acyl = subs[0].0.clone();
            let alkyl = subs[1].0.clone();
            Some((acyl, alkyl))
        } else {
            None
        }
    }

    /// Match an ether functional group in the given graph.
    ///
    /// This function searches for an oxygen atom that is connected via single bonds
    /// to exactly two carbon atoms. When found, the oxygen is removed (along with its edges)
    /// and the remaining graph is partitioned into two connected fragments.
    fn match_ether(graph: &MoleculeGraph) -> Option<(MoleculeGraph, MoleculeGraph)> {
        let (subs, _) = disconnect_subgraph_by_pattern(graph, &Substituent::ether()).ok()?;
        if subs.len() == 2 {
            let alkyl1 = subs[0].0.clone();
            let alkyl2 = subs[1].0.clone();
            Some((alkyl1, alkyl2))
        } else {
            None
        }

        // for node in graph.node_indices() {
        //     if graph[node].kind == ElementType::O {
        //         let neighbors: Vec<_> = graph.neighbors(node).collect();
        //         if neighbors.len() == 2 &&
        //            neighbors.iter().all(|&n| graph[n].kind == ElementType::C)
        //         {
        //             // Remove the oxygen node from the graph to disconnect the two sides.
        //             let g_without_o = remove_node_from_graph(graph, node);
        //             // Compute the connected components of the resulting graph.
        //             let comps = connected_components(&g_without_o);
        //             if comps.len() == 2 {
        //                 let g1 = subgraph_from_component(&g_without_o, &comps[0]);
        //                 let g2 = subgraph_from_component(&g_without_o, &comps[1]);
        //                 return Some((g1, g2));
        //             }
        //         }
        //     }
        // }
        // None
    }

    /// Match a carboxylic acid pattern in the graph.
    ///
    /// We look for a carbon that is doubly bonded to an oxygen and singly bonded to an oxygen
    /// that is terminal (i.e. has only one neighbor). The terminal oxygen is then removed by
    /// splitting the graph along that bond, and the fragment not containing the carbon is treated
    /// as the target group (for example, H for formic acid).
    fn match_carboxylic_acid(graph: &MoleculeGraph) -> Option<MoleculeGraph> {
        let (subs, _) = disconnect_subgraph_by_pattern(graph, &Substituent::carboxyl_group()).ok()?;
        if subs.len() == 1 {
            let target = subs[0].0.clone();
            Some(target)
        } else {
            None
        }
    }

    pub fn to_molecule_graph(&self) -> MoleculeGraph {
        let mut clone = self.clone();
        clone.merge_all_substituents().to_molecule_graph()
    }

    fn count_r_groups(&self) -> usize {
        match self {
            OrganicMolecule::Substituent(substituent) => substituent.count_r_groups(),
            OrganicMolecule::Compound(graph) => {
                let total_r_groups = graph.node_indices().map(|idx| graph[idx].count_r_groups().saturating_sub(graph.neighbors(idx).count())).sum::<usize>();
                // total_r_groups - graph.edge_count()
                total_r_groups
            }
        }
    }

    pub fn visualize(&self, filename: &str) -> Result<()> {
        let mut graph = self.to_molecule_graph();
        graph.hydrogenate();
        visualize_graph(&graph, "tmp.dot", Some(filename)).map_err(|e| anyhow!("Failed to visualize: {}", e))?;
        Ok(())
    }

    pub fn alkane(n: usize, branches: impl IntoIterator<Item = (usize, Self)>) -> Self {
        let branches = branches.into_iter().collect::<Vec<_>>();
        let mut graph = UnGraph::new_undirected();

        // Sort the branches by position
        let mut branches = branches;
        branches.sort_by_key(|(pos, _)| *pos);

        // Create the main chain
        let main_chain = Substituent::alkyl_group(n, branches.iter().map(|(pos, _)| *pos));
        let main_chain_node = graph.add_node(Self::from(main_chain));

        // Add the branches
        for (pos, branch) in branches {
            let branch_node = graph.add_node(branch);
            graph.add_edge(main_chain_node, branch_node, pos);
        }

        OrganicMolecule::Compound(graph).check_r_groups(0..=1, "Resulting alkane must have 0 or 1 R groups")
    }

    pub fn alkyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        let alkyl = Substituent::alkyl_group(n, r_group_positions);
        OrganicMolecule::Substituent(alkyl).check_r_groups(0..=1, "Resulting alkyl group must have 0 or 1 R groups")
    }

    pub fn methane() -> Self {
        Self::alkane(1, vec![])
    }

    pub fn ethane() -> Self {
        Self::alkane(2, vec![])
    }

    pub fn propane() -> Self {
        Self::alkane(3, vec![])
    }

    pub fn butane() -> Self {
        Self::alkane(4, vec![])
    }

    pub fn pentane() -> Self {
        Self::alkane(5, vec![])
    }

    pub fn hexane() -> Self {
        Self::alkane(6, vec![])
    }

    pub fn methyl() -> Self {
        Self::alkyl_group(1, vec![1])
    }

    pub fn ethyl() -> Self {
        Self::alkyl_group(2, vec![1])
    }

    pub fn propyl() -> Self {
        Self::alkyl_group(3, vec![1])
    }

    pub fn butyl() -> Self {
        Self::alkyl_group(4, vec![1])
    }

    pub fn pentyl() -> Self {
        Self::alkyl_group(5, vec![1])
    }

    pub fn hexyl() -> Self {
        Self::alkyl_group(6, vec![1])
    }
    
    pub fn carboxylic_acid(target: impl Into<Self>) -> Self {
        let mut graph = UnGraph::new_undirected();
        
        let target = target.into();
        let carboxylic_acid = Substituent::carboxyl_group().into();
        let carboxylic_acid_node = graph.add_node(carboxylic_acid);
        let target_node = graph.add_node(target);
        graph.add_edge(carboxylic_acid_node, target_node, 0);

        OrganicMolecule::Compound(graph).check_r_groups(0..=1, "Resulting carboxylic acid must have 0 or 1 R groups")
    }

    pub fn cyclic_alkane(n: usize, branches: impl IntoIterator<Item = (usize, Self)>) -> Self {
        let branches = branches.into_iter().collect::<Vec<_>>();
        let mut graph = UnGraph::new_undirected();

        // Sort the branches by position
        let mut branches = branches;
        branches.sort_by_key(|(pos, _)| *pos);

        // Create the main chain
        let main_chain = Substituent::cyclic_alkyl_group(n, branches.iter().map(|(pos, _)| *pos));
        let main_chain_node = graph.add_node(Self::from(main_chain));

        // Add the branches
        for (pos, branch) in branches {
            let branch_node = graph.add_node(branch);
            graph.add_edge(main_chain_node, branch_node, pos);
        }

        OrganicMolecule::Compound(graph).check_r_groups(0..=1, "Resulting cyclic alkane must have 0 or 1 R groups")
    }

    fn add_r_group(&mut self) {
        if let OrganicMolecule::Substituent(sub) = self {
            warn!("Adding R group to substituent with no R groups: {:?}", sub);
            // Connect an R group to the first carbon
            let carbon = sub.0.node_indices().find(|&node| sub.0[node].is_carbon()).unwrap();
            let r_group = sub.0.add_node(Element::r_group(0));
            sub.0.add_edge(carbon, r_group, Bond::Single);
        } else {
            // Topologically sort the graph
            let graph = self.mut_graph();
            let sorted = petgraph::algo::toposort(&*graph, None).unwrap();
            // Get the first node
            let first_node = sorted.first().unwrap();
            // Add an R group to the first node
            graph[*first_node].add_r_group();
        }
    }

    pub fn cyclic_alkyl_group(n: usize, r_group_positions: impl IntoIterator<Item = usize>) -> Self {
        let alkyl = Substituent::cyclic_alkyl_group(n, r_group_positions);
        OrganicMolecule::Substituent(alkyl).check_r_groups(0..=1, "Resulting cyclic alkyl group must have 0 or 1 R groups")
    }

    pub fn phenyl() -> Self {
        Substituent::phenyl_group().into()
    }

    pub fn halide(elem: ElementType) -> Self {
        Substituent::halide(elem).into()
    }

    pub fn hydrogen() -> Self {
        Substituent::hydrogen().into()
    }

    pub fn ether(alkyl1: impl Into<Self>, alkyl2: impl Into<Self>) -> Self {
        let mut graph = UnGraph::new_undirected();
        // Add the ether
        let ether = graph.add_node(Substituent::ether().into());

        let alkyl1 = graph.add_node(alkyl1.into());
        let alkyl2 = graph.add_node(alkyl2.into());

        graph.add_edge(ether, alkyl1, 0);
        graph.add_edge(ether, alkyl2, 1);

        OrganicMolecule::Compound(graph).check_r_groups(0..=2, "Resulting ether must have 0, 1, or 2 R groups")
    }

    pub fn ester(acyl: impl Into<Self>, alkyl: impl Into<Self>) -> Self {
        let mut graph = UnGraph::new_undirected();
        // Add the ether
        let ester = graph.add_node(Substituent::ester().into());

        let acyl = graph.add_node(acyl.into());
        let alkyl = graph.add_node(alkyl.into());

        graph.add_edge(ester, acyl, 0);
        graph.add_edge(ester, alkyl, 1);

        OrganicMolecule::Compound(graph).check_r_groups(0..=2, "Resulting ether must have 0, 1, or 2 R groups")
    }

    pub fn alcohol(alkyl: impl Into<Self>) -> Self {
        let mut graph = UnGraph::new_undirected();
        let alkyl = alkyl.into();
        graph.add_node(alkyl);
        OrganicMolecule::Compound(graph).check_r_groups(0..=1, "Alcohol must have 0 or 1 R groups")
    }

    fn check_r_groups(self, expected: impl RangeBounds<usize> + Debug, msg: impl ToString) -> Self {
        let r_groups = self.count_r_groups();
        if !expected.contains(&r_groups) {
            panic!("{}: Expected {:?} R groups, found {}", msg.to_string(), expected, r_groups);
        }
        self
    }

    fn merge_substituents(&mut self, s1_node_idx: NodeIndex, s2_node_idx: NodeIndex) {
        if let OrganicMolecule::Compound(ref mut graph) = self {
            let mut s1_orgo = graph[s1_node_idx].clone();
            let mut s2_orgo = graph[s2_node_idx].clone();
            
            if s1_orgo.count_r_groups() < 1 {
                s1_orgo.add_r_group();
            }
            if s2_orgo.count_r_groups() < 1 {
                s2_orgo.add_r_group();
            }

            let mut s1 = s1_orgo.as_substituent().clone();
            let s2 = s2_orgo.as_substituent().clone();

            // Connect the two substituents
            s1.connect(&s2).expect("Failed to connect substituents");

            // Store the new substituent back in the graph
            graph[s1_node_idx] = s1.into();
    
            // Remove the other substituent
            graph.remove_node(s2_node_idx);
        }
    }

    fn as_substituent(&self) -> Substituent {
        match self {
            OrganicMolecule::Substituent(substituent) => substituent.clone(),
            Self::Compound(_) => {
                self.clone().merge_all_substituents()
            }
        }
    }

    fn is_compound(&self) -> bool {
        match self {
            OrganicMolecule::Compound(_) => true,
            _ => false
        }
    }

    fn substituent_count(&self) -> usize {
        match self {
            OrganicMolecule::Compound(graph) => graph.node_count(),
            _ => 1
        }
    }

    fn ref_graph(&self) -> &UnGraph<Self, usize> {
        match self {
            OrganicMolecule::Compound(graph) => graph,
            _ => unreachable!()
        }
    }

    fn mut_graph(&mut self) -> &mut UnGraph<Self, usize> {
        match self {
            OrganicMolecule::Compound(graph) => graph,
            _ => unreachable!()
        }
    }

    fn merge_all_substituents(&mut self) -> Substituent {
        if self.is_compound() {
            while self.substituent_count() > 1 {
                // Sort the edges by weight
                let mut edges: Vec<_> = self.ref_graph().edge_indices().collect();
                edges.sort_by_key(|&edge| self.ref_graph()[edge]);
                
                // Get the two substituents to connect
                let (s1_node_idx, s2_node_idx) = self.ref_graph().edge_endpoints(edges[0]).unwrap();
                
                // Merge the two substituents
                self.merge_substituents(s1_node_idx, s2_node_idx);
            }
            self.ref_graph()[self.ref_graph().node_indices().next().unwrap()].as_substituent()
        } else {
            self.as_substituent().clone()
        }
    }

    pub fn to_smiles(&self) -> Result<String> {
        Ok(molecule_to_smiles(&self.to_molecule_graph()))
    }

    pub fn to_iupac(&self) -> Result<String> {
        Ok(iupac_name(&self.to_molecule_graph()))
    }
}


impl Debug for OrganicMolecule {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            OrganicMolecule::Substituent(substituent) => write!(f, "{:?}", substituent),
            OrganicMolecule::Compound(graph) => {
                let mut fmt = f.debug_map();
                for node in graph.node_indices() {
                    fmt.entry(&node, &graph[node]);
                }
                fmt.finish()
            }
        }
    }
}

impl From<Substituent> for OrganicMolecule {
    fn from(substituent: Substituent) -> Self {
        OrganicMolecule::Substituent(substituent)
    }
}

impl From<Element> for OrganicMolecule {
    fn from(elem: Element) -> Self {
        OrganicMolecule::from(Substituent::from(elem))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_smiles(molecule: &OrganicMolecule) {
        // Convert the molecule to a SMILES string
        let smiles = molecule.to_smiles().unwrap();

        // Parse the SMILES string back into a molecule
        let parsed = OrganicMolecule::parse_smiles(&smiles).unwrap();
        parsed.visualize("parsed.png").unwrap();
        // Check that the parsed molecule is the same as the original
        assert!(molecule.is_same_as(&parsed), "Molecule mismatch: {:?} != {:?}", molecule, parsed);
    }

    #[test]
    fn test_smiles() {
        init_logging("info");
        // let hexane = OrganicMolecule::parse_smiles("CCCCCC").unwrap();
        // println!("{:#?}", hexane);

        // let cyclohexane = OrganicMolecule::parse_smiles("C1CCCCC1").unwrap();
        // println!("{:#?}", cyclohexane);

        // let methylcyclopropane = OrganicMolecule::parse_smiles("CC1CC1").unwrap();
        // println!("{:#?}", methylcyclopropane);

        // let _4ethyl_2methylhexane = OrganicMolecule::parse_smiles("CCC(C)C(C)CC").unwrap();
        // println!("{:#?}", _4ethyl_2methylhexane);
        // _4ethyl_2methylhexane.visualize("4ethyl-2methylhexane.png").unwrap();

        info!("benzyl acetate");
        let benzyl_acetate = OrganicMolecule::parse_smiles("c1ccccc1COC(=O)C").unwrap();
        benzyl_acetate.visualize("benzyl_acetate.png").unwrap();
        info!("benzyl acetate {:?}", benzyl_acetate);

        info!("isoamyl acetate");
        let isoamyl_acetate = OrganicMolecule::parse_smiles("O=C(OCCC(C)C)C").unwrap();
        isoamyl_acetate.visualize("isoamyl_acetate.png").unwrap();
        info!("isoamyl acetate: {:?}", isoamyl_acetate);

        info!("ethyl acetate");
        let isobutyl_formate = OrganicMolecule::parse_smiles("O=COCC(C)C").unwrap();
        isobutyl_formate.visualize("isobutyl_formate.png").unwrap();
        info!("ethyl acetate: {:?}", isobutyl_formate);

        info!("ethyl benzoate");
        let ethyl_benzoate = OrganicMolecule::parse_smiles("O=C(OCC)c1ccccc1").unwrap();
        ethyl_benzoate.visualize("ethyl_benzoate.png").unwrap();
        info!("ethyl benzoate: {:?}", ethyl_benzoate);

        // check_smiles(&benzyl_acetate);
        // check_smiles(&isoamyl_acetate);
        // check_smiles(&isobutyl_formate);
        check_smiles(&ethyl_benzoate);
    }

    #[test]
    fn test_organic_molecule() -> Result<()> {
        let hexane = OrganicMolecule::alkane(6, vec![]);
        println!("{:#?}", hexane);
        hexane.visualize("hexane.png")?;

        let cyclohexane = OrganicMolecule::cyclic_alkane(6, vec![]);
        println!("{:#?}", cyclohexane);
        cyclohexane.visualize("cyclohexane.png")?;


        let methylcyclopropane = OrganicMolecule::cyclic_alkane(3,
            vec![(1, OrganicMolecule::methyl())]);
        println!("{:#?}", methylcyclopropane);
        methylcyclopropane.visualize("methylcyclopropane.png")?;

        let _4ethyl_2methylhexane = OrganicMolecule::alkane(6, [
                (2, OrganicMolecule::methyl()),
                (4, OrganicMolecule::ethyl()),
            ]);
        println!("{:#?}", _4ethyl_2methylhexane);
        _4ethyl_2methylhexane.visualize("4ethyl-2methylhexane.png")?;

        let _2_3_5trimethyl_4_propylheptane = OrganicMolecule::alkane(7, [
            (4, OrganicMolecule::propyl()),
            (2, OrganicMolecule::methyl()),
            (3, OrganicMolecule::methyl()),
            (5, OrganicMolecule::methyl()),
        ]);
        println!("{:#?}", _2_3_5trimethyl_4_propylheptane);
        _2_3_5trimethyl_4_propylheptane.visualize("2,3,5trimethyl-4propylheptane.png")?;
    
        let formic_acid = OrganicMolecule::carboxylic_acid(Element::H);
        println!("{:#?}", formic_acid);
        formic_acid.visualize("formic_acid.png")?;

        let acetic_acid = OrganicMolecule::carboxylic_acid(OrganicMolecule::methyl());
        println!("{:#?}", acetic_acid);
        
        Ok(())
    }

    #[test]
    fn test_substituents() -> Result<()> {
        let hexane = Substituent::alkyl_group(6, vec![1, 3, 5]);
        println!("{:#?}", hexane);

        hexane.visualize("hexane.png")?;

        let methylpentane = Substituent::alkyl_group(5, vec![3, 3])
            .then_connect(Substituent::alkyl_group(1, vec![1]))?
            .then_connect(Substituent::alkyl_group(1, vec![1]))?;
    
        println!("{:#?}", methylpentane);

        methylpentane.visualize("3-methylpentane.png")?;

        let cyclohexane = Substituent::cyclic_alkyl_group(6, vec![1, 3, 5])
            .then_connect(Substituent::alkyl_group(3, vec![1]))?
            .then_connect(Substituent::alkyl_group(3, vec![1]))?
            .then_connect(Substituent::alkyl_group(3, vec![1]))?;
        println!("{:#?}", cyclohexane);
        cyclohexane.visualize("cyclohexane.png")?;

        let isobutyl = Substituent::alkyl_group(3, vec![2])
            .then_connect(Substituent::alkyl_group(1, vec![1, 1]))?;

        isobutyl.visualize("isobutyl.png")?;

        let sec_butyl = Substituent::alkyl_group(3, vec![1, 2])
            .then_connect(Substituent::alkyl_group(1, vec![1]))?;
        sec_butyl.visualize("sec_butyl.png")?;

        let tert_butyl = Substituent::alkyl_group(1, vec![1, 1, 1, 1])
            .then_connect(Substituent::alkyl_group(1, vec![1]))?
            .then_connect(Substituent::alkyl_group(1, vec![1]))?
            .then_connect(Substituent::alkyl_group(1, vec![1]))?;
        tert_butyl.visualize("tert_butyl.png")?;

        let carboxyl = Substituent::carboxyl_group();
        carboxyl.visualize("carboxyl.png")?;

        let methylpropanoate = Substituent::ester()
            .then_connect(Substituent::alkyl_group(2, vec![1]))?
            .then_connect(Substituent::alkyl_group(1, vec![1]))?;
        methylpropanoate.visualize("methylpropanoate.png")?;

        let chloro = Substituent::halide(ElementType::Cl);
        chloro.visualize("cloro.png")?;

        let _3cloro_2methylpentane = Substituent::alkyl_group(5, vec![2, 3])
            .then_connect(Substituent::alkyl_group(1, vec![1]))?
            .then_connect(Substituent::halide(ElementType::Cl))?;

        _3cloro_2methylpentane.visualize_raw("3cloro-2methylpentane.png")?;

        let _2bromo_3methylbutane = Substituent::alkyl_group(4, vec![2, 3])
            .then_connect(Substituent::alkyl_group(1, vec![1]))?
            .then_connect(Substituent::halide(ElementType::Br))?;
        _2bromo_3methylbutane.visualize_raw("2bromo-3methylbutane.png")?;

        Ok(())
    }
}