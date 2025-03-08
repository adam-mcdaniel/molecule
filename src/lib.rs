use anyhow::{Context, Result};
use parse::SmilesError;
use petgraph::dot::{Config, Dot};
use petgraph::graph::EdgeReference;
use petgraph::graph::NodeIndex;
use petgraph::graph::{EdgeIndex, Graph, UnGraph};
use petgraph::visit::EdgeRef;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs::File;
use std::io::Write;
use thiserror::Error as ThisError;
use tracing::*;

#[derive(ThisError, Debug)]
pub enum ParseError {
    #[error("Unknown element: {0}")]
    UnknownElement(String),

    #[error("Invalid SMILES string: {0}")]
    SmilesError(SmilesError),
}

#[derive(ThisError, Debug)]
pub enum Error {
    #[error("Parse error: {0}")]
    ParseError(ParseError),
}

impl From<ParseError> for Error {
    fn from(e: ParseError) -> Self {
        Error::ParseError(e)
    }
}

mod parse;
pub use parse::*;

mod visualize;
pub use visualize::*;

mod naming;
pub use naming::*;

mod organic;
pub use organic::*;

pub fn init_logging(level: &str) {
    let _ = tracing_subscriber::fmt()
        .with_max_level(match level {
            "trace" => tracing::metadata::LevelFilter::TRACE,
            "debug" => tracing::metadata::LevelFilter::DEBUG,
            "info" => tracing::metadata::LevelFilter::INFO,
            "warn" => tracing::metadata::LevelFilter::WARN,
            "error" => tracing::metadata::LevelFilter::ERROR,
            _ => tracing::metadata::LevelFilter::OFF,
        })
        .without_time()
        // .with_target(false)
        .try_init();
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ElementType {
    C,
    H,
    O,
    N,
    F,
    Cl,
    Br,
    S,
    I,
    RGroup(usize),
}

impl ElementType {
    pub fn as_aromatic(self) -> Element {
        Element::aromatic(self)
    }

    pub fn as_aliphatic(self) -> Element {
        Element::aliphatic(self)
    }

    pub fn available_single_bonds(&self) -> Option<usize> {
        Some(match self {
            ElementType::C => 4,
            ElementType::H => 1,
            ElementType::O => 2,
            ElementType::N => 3,
            ElementType::F => 1,
            ElementType::Cl => 1,
            ElementType::Br => 1,
            ElementType::S => 2,
            ElementType::I => 1,
            ElementType::RGroup(_) => return None,
        })
    }

    pub fn diameter(&self) -> f64 {
        match self {
            ElementType::C => 1.54,
            ElementType::H => 1.2,
            ElementType::O => 1.52,
            ElementType::N => 1.47,
            ElementType::F => 1.47,
            ElementType::Cl => 1.75,
            ElementType::Br => 1.85,
            ElementType::S => 1.8,
            ElementType::I => 1.98,
            ElementType::RGroup(_) => 1.54,
        }
    }

    pub fn is_carbon(&self) -> bool {
        matches!(self, ElementType::C)
    }

    pub fn is_hydrogen(&self) -> bool {
        matches!(self, ElementType::H)
    }

    pub fn is_oxygen(&self) -> bool {
        matches!(self, ElementType::O)
    }

    pub fn is_nitrogen(&self) -> bool {
        matches!(self, ElementType::N)
    }

    pub fn is_sulfur(&self) -> bool {
        matches!(self, ElementType::S)
    }

    pub fn is_iodine(&self) -> bool {
        matches!(self, ElementType::I)
    }

    pub fn is_r_group(&self) -> bool {
        matches!(self, ElementType::RGroup(_))
    }

    pub fn symbol(&self) -> String {
        match self {
            ElementType::C => "C",
            ElementType::H => "H",
            ElementType::O => "O",
            ElementType::N => "N",
            ElementType::F => "F",
            ElementType::Cl => "Cl",
            ElementType::Br => "Br",
            ElementType::S => "S",
            ElementType::I => "I",
            ElementType::RGroup(count) => return "R".to_string() + "'".repeat(*count).as_str(),
        }
        .to_string()
    }

    pub fn is_halogen(&self) -> bool {
        matches!(
            self,
            ElementType::F | ElementType::Cl | ElementType::Br | ElementType::I
        )
    }

    pub fn is_heteroatom(&self) -> bool {
        !self.is_carbon()
    }
}

impl From<Element> for ElementType {
    fn from(element: Element) -> Self {
        element.kind
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    pub kind: ElementType,
    pub aromatic: bool,
}

impl Element {
    pub const H: Self = Self {
        kind: ElementType::H,
        aromatic: false,
    };

    pub const C: Self = Self {
        kind: ElementType::C,
        aromatic: false,
    };

    pub const O: Self = Self {
        kind: ElementType::O,
        aromatic: false,
    };

    pub const N: Self = Self {
        kind: ElementType::N,
        aromatic: false,
    };

    pub const F: Self = Self {
        kind: ElementType::F,
        aromatic: false,
    };

    #[allow(non_upper_case_globals)]
    pub const Cl: Self = Self {
        kind: ElementType::Cl,
        aromatic: false,
    };

    #[allow(non_upper_case_globals)]
    pub const Br: Self = Self {
        kind: ElementType::Br,
        aromatic: false,
    };

    pub const S: Self = Self {
        kind: ElementType::S,
        aromatic: false,
    };

    pub const I: Self = Self {
        kind: ElementType::I,
        aromatic: false,
    };

    pub fn available_single_bonds(&self) -> Option<usize> {
        self.kind.available_single_bonds()
    }

    pub fn r_group(n: usize) -> Self {
        Self {
            kind: ElementType::RGroup(n),
            aromatic: false,
        }
    }

    pub fn as_aromatic(self) -> Self {
        Self {
            aromatic: true,
            ..self
        }
    }

    pub fn as_aliphatic(self) -> Self {
        Self {
            aromatic: false,
            ..self
        }
    }

    pub fn from_smiles(smiles: impl AsRef<str>) -> anyhow::Result<Self> {
        Ok(match smiles.as_ref() {
            "H" => Self::H,
            "C" => Self::C,
            "O" => Self::O,
            "N" => Self::N,
            "F" => Self::F,
            "Cl" => Self::Cl,
            "Br" => Self::Br,
            "S" => Self::S,

            "c" => Self::C.as_aromatic(),
            "o" => Self::O.as_aromatic(),
            "n" => Self::N.as_aromatic(),
            "s" => Self::S.as_aromatic(),
            // Parse R groups
            // The number of apostrophes indicates the Nth R group
            "R" => Self::r_group(1),
            "R'" => Self::r_group(2),
            "R''" => Self::r_group(3),
            "R'''" => Self::r_group(4),

            unknown => {
                return Err(ParseError::UnknownElement(unknown.to_string()))
                    .context("While parsing SMILES string")
            }
        })
    }

    pub fn aromatic(kind: ElementType) -> Self {
        Self {
            kind,
            aromatic: true,
        }
    }

    pub fn aliphatic(kind: ElementType) -> Self {
        Self {
            kind,
            aromatic: false,
        }
    }

    pub fn is_r_group(&self) -> bool {
        self.kind.is_r_group()
    }

    pub fn is_carbon(&self) -> bool {
        self.kind.is_carbon()
    }

    pub fn is_hydrogen(&self) -> bool {
        self.kind.is_hydrogen()
    }

    pub fn is_oxygen(&self) -> bool {
        self.kind.is_oxygen()
    }

    pub fn is_nitrogen(&self) -> bool {
        self.kind.is_nitrogen()
    }

    pub fn is_sulfur(&self) -> bool {
        self.kind.is_sulfur()
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic
    }

    pub fn is_aliphatic(&self) -> bool {
        !self.aromatic
    }

    pub fn is_halogen(&self) -> bool {
        self.kind.is_halogen()
    }

    pub fn is_heteroatom(&self) -> bool {
        self.kind.is_heteroatom()
    }

    pub fn symbol(&self) -> String {
        self.kind.symbol()
    }

    pub fn diameter(&self) -> f64 {
        self.kind.diameter()
    }
}

impl PartialEq<ElementType> for Element {
    fn eq(&self, other: &ElementType) -> bool {
        self.kind == *other
    }
}

impl From<ElementType> for Element {
    fn from(kind: ElementType) -> Self {
        Self {
            kind,
            aromatic: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Bond {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl Bond {
    pub fn order(&self) -> f64 {
        match self {
            Bond::Single => 1.0,
            Bond::Double => 2.0,
            Bond::Triple => 3.0,
            Bond::Aromatic => 1.5,
        }
    }

    pub fn is_aromatic(&self) -> bool {
        matches!(self, Bond::Aromatic)
    }
}

pub type MoleculeGraph = petgraph::graph::UnGraph<Element, Bond>;

impl From<Element> for MoleculeGraph {
    fn from(element: Element) -> Self {
        let mut graph = MoleculeGraph::default();
        graph.add_node(element);
        graph
    }
}

pub(crate) fn connect(
    dst: &mut MoleculeGraph,
    src: &MoleculeGraph,
    dst_node: NodeIndex,
    src_node: NodeIndex,
    bond: Bond,
) {
    // Take all the nodes and edges from the source graph
    let mut node_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    for node in src.node_indices() {
        let new_node = dst.add_node(src[node]);
        node_map.insert(node, new_node);
    }

    for edge in src.edge_indices() {
        let (source, target) = src.edge_endpoints(edge).unwrap();

        let source = node_map[&source];
        let target = node_map[&target];
        dst.add_edge(source, target, *src.edge_weight(edge).unwrap());
    }

    // Add the connection between the two graphs
    let new_src_node = node_map[&src_node];
    dst.add_edge(dst_node, new_src_node, bond);
}

fn hydrogenate(graph: &mut MoleculeGraph) {
    for node in graph.node_indices() {
        let element = graph[node];
        if let Some(available_bonds) = element.available_single_bonds() {
            let mut available_bonds = available_bonds as f64;

            for edge in graph.edges(node) {
                available_bonds -= edge.weight().order();
            }

            if available_bonds < 0.0 {
                warn!(
                    "Found element {} with negative available bonds in {graph:#?}",
                    element.symbol()
                );
            }

            for _ in 0..available_bonds.ceil() as usize {
                let h_node = graph.add_node(Element::H);
                graph.add_edge(node, h_node, Bond::Single);
            }
        }
    }
}

pub trait Hydrogenate {
    fn hydrogenate(&mut self);
}

impl Hydrogenate for MoleculeGraph {
    fn hydrogenate(&mut self) {
        hydrogenate(self);
    }
}

/// Given a MoleculeGraph, returns its molecular formula as a string (e.g. "C9H10O2").
/// For organic formulas, we print carbon (C) first and hydrogen (H) second (if present),
/// and then all other elements in alphabetical order.
pub fn formula_string(graph: &MoleculeGraph) -> String {
    // Count each element symbol.
    let mut counts: HashMap<String, usize> = HashMap::new();
    for node in graph.node_indices() {
        // Use the Element::symbol() method.
        let symbol = graph[node].symbol();
        *counts.entry(symbol).or_insert(0) += 1;
    }

    let mut output = String::new();
    // Print Carbon first.
    if let Some(&c_count) = counts.get("C") {
        output.push_str("C");
        if c_count > 1 {
            write!(output, "{}", c_count).unwrap();
        }
        counts.remove("C");
    }
    // Then Hydrogen.
    if let Some(&h_count) = counts.get("H") {
        output.push_str("H");
        if h_count > 1 {
            write!(output, "{}", h_count).unwrap();
        }
        counts.remove("H");
    }
    // Print the remaining elements in alphabetical order.
    for (symbol, count) in counts.iter() {
        output.push_str(symbol);
        if *count > 1 {
            write!(output, "{}", count).unwrap();
        }
    }
    output
}

/// A helper structure to record a ring closure edge.
#[derive(Debug, Clone)]
struct RingClosure {
    opening: NodeIndex, // the ancestor where the ring is opened
    closing: NodeIndex, // the descendant where the ring is closed
    bond: Bond,         // the bond used for the closure
    digit: usize,       // assigned ring closure digit
}

/// First pass: Compute a spanning tree of the graph (via DFS) and record ring closures.
/// Returns a mapping from each node (except the root) to its parent in the DFS tree,
/// and a vector of RingClosure records for every back edge encountered.
fn compute_spanning_tree_and_ring_closures(
    graph: &MoleculeGraph,
    start: NodeIndex,
) -> (HashMap<NodeIndex, NodeIndex>, Vec<RingClosure>) {
    let mut parent: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    let mut ring_closures: Vec<RingClosure> = Vec::new();
    let mut visited: HashSet<NodeIndex> = HashSet::new();
    let mut path: Vec<NodeIndex> = Vec::new(); // current DFS stack

    fn dfs(
        graph: &MoleculeGraph,
        current: NodeIndex,
        parent: &mut HashMap<NodeIndex, NodeIndex>,
        ring_closures: &mut Vec<RingClosure>,
        visited: &mut HashSet<NodeIndex>,
        path: &mut Vec<NodeIndex>,
    ) {
        visited.insert(current);
        path.push(current);
        // Iterate over neighbors (in sorted order for stability).
        let mut nbrs: Vec<NodeIndex> = graph.neighbors(current).collect();
        nbrs.sort();
        for nbr in nbrs {
            // Skip the parent edge (if any)
            if let Some(&p) = parent.get(&current) {
                if nbr == p {
                    continue;
                }
            }
            if !visited.contains(&nbr) {
                parent.insert(nbr, current);
                dfs(graph, nbr, parent, ring_closures, visited, path);
            } else if path.contains(&nbr) {
                // Found a back edge: current -> nbr, where nbr is an ancestor.
                // Record it as a ring closure. We assign a digit later.
                // Only record if the edge isn’t already recorded (using ordered key).
                let key = if current < nbr {
                    (current, nbr)
                } else {
                    (nbr, current)
                };
                // Check if already recorded:
                if !ring_closures.iter().any(|rc| {
                    let existing_key = if rc.opening < rc.closing {
                        (rc.opening, rc.closing)
                    } else {
                        (rc.closing, rc.opening)
                    };
                    existing_key == key
                }) {
                    ring_closures.push(RingClosure {
                        opening: nbr,
                        closing: current,
                        bond: *graph
                            .edge_weight(graph.find_edge(current, nbr).unwrap())
                            .unwrap(),
                        digit: 0, // placeholder; will assign digits in a later pass
                    });
                }
            }
        }
        path.pop();
    }

    dfs(
        graph,
        start,
        &mut parent,
        &mut ring_closures,
        &mut visited,
        &mut path,
    );
    (parent, ring_closures)
}
/// Second pass: Generate SMILES from the spanning tree, inserting ring–closure markers
/// and forcing branch notation when the current atom is a branching center.
fn generate_smiles_from_tree(
    graph: &MoleculeGraph,
    root: NodeIndex,
    parent: &HashMap<NodeIndex, NodeIndex>,
    ring_closures: &mut Vec<RingClosure>,
) -> String {
    // First, assign digits to ring closures.
    ring_closures.sort_by_key(|rc| (rc.opening.index(), rc.closing.index()));
    for (i, rc) in ring_closures.iter_mut().enumerate() {
        rc.digit = i + 1;
    }

    // Build the children mapping from the parent mapping.
    let mut children: HashMap<NodeIndex, Vec<NodeIndex>> = HashMap::new();
    for (&child, &p) in parent.iter() {
        children.entry(p).or_default().push(child);
    }
    // Sort children for stable output.
    for child_list in children.values_mut() {
        child_list.sort();
    }

    // Build lookups for ring closures.
    let mut open_map: HashMap<NodeIndex, Vec<&RingClosure>> = HashMap::new();
    let mut close_map: HashMap<NodeIndex, Vec<&RingClosure>> = HashMap::new();
    for rc in ring_closures.iter() {
        open_map.entry(rc.opening).or_default().push(rc);
        close_map.entry(rc.closing).or_default().push(rc);
    }

    // Recursively generate SMILES with a modified DFS.
    fn dfs_tree(
        graph: &MoleculeGraph,
        current: NodeIndex,
        children: &HashMap<NodeIndex, Vec<NodeIndex>>,
        open_map: &HashMap<NodeIndex, Vec<&RingClosure>>,
        close_map: &HashMap<NodeIndex, Vec<&RingClosure>>,
        is_root: bool,
    ) -> String {
        // Start with the current atom's symbol.
        let mut s = graph[current].symbol();

        // Insert ring closure markers that open at this atom.
        if let Some(rcs) = open_map.get(&current) {
            for rc in rcs {
                let bond = match rc.bond {
                    Bond::Single => "",
                    Bond::Double => "=",
                    Bond::Triple => "#",
                    Bond::Aromatic => ":",
                };
                s.push_str(bond);
                s.push_str(&format_ring(rc.digit));
            }
        }
        // Insert ring closures that close at this atom after the atom symbol.
        if let Some(rcs) = close_map.get(&current) {
            for rc in rcs {
                let bond = match rc.bond {
                    Bond::Single => "",
                    Bond::Double => "=",
                    Bond::Triple => "#",
                    Bond::Aromatic => ":",
                };
                s.push_str(bond);
                s.push_str(&format_ring(rc.digit));
            }
        }
        // Determine how many children this node has in the spanning tree.
        let child_list = children.get(&current);
        let num_children = child_list.map_or(0, |v| v.len());
        // Compute the total degree in the spanning tree for a non-root node.
        // For the root, we force branch notation for every child.
        let force_branch = if is_root { true } else { num_children > 1 };

        // Process children.
        if let Some(child_nodes) = child_list {
            let mut branch_str = String::new();
            if force_branch {
                // For each child, output as a branch (with bond symbol) enclosed in parentheses.
                for &child in child_nodes.iter() {
                    branch_str.push_str(&format!(
                        "({}{})",
                        bond_str(graph, current, child),
                        dfs_tree(graph, child, children, open_map, close_map, false)
                    ));
                }
            } else {
                // Only one child: inline the continuation.
                let child = child_nodes[0];
                branch_str.push_str(&bond_str(graph, current, child));
                branch_str.push_str(&dfs_tree(
                    graph, child, children, open_map, close_map, false,
                ));
            }
            s.push_str(&branch_str);
        }
        s
    }

    dfs_tree(graph, root, &children, &open_map, &close_map, true)
}

/// Returns the bond symbol for the edge between nodes `a` and `b`.
fn bond_str(graph: &MoleculeGraph, a: NodeIndex, b: NodeIndex) -> String {
    if let Some(edge) = graph.find_edge(a, b) {
        let bond = graph.edge_weight(edge).unwrap();
        match bond {
            Bond::Single => "".to_string(),
            Bond::Double => "=".to_string(),
            Bond::Triple => "#".to_string(),
            Bond::Aromatic => ":".to_string(),
        }
    } else {
        "".to_string()
    }
}

/// Formats a ring closure digit according to SMILES rules.
fn format_ring(digit: usize) -> String {
    if digit < 10 {
        digit.to_string()
    } else {
        format!("%{}", digit)
    }
}

/// Top-level function to convert a MoleculeGraph into a SMILES string.
pub fn molecule_to_smiles(graph: &MoleculeGraph) -> String {
    let start = graph.node_indices().next().unwrap();
    let (parent, mut ring_closures) = compute_spanning_tree_and_ring_closures(graph, start);
    generate_smiles_from_tree(graph, start, &parent, &mut ring_closures)
}
