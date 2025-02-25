
use petgraph::dot::{Config, Dot};
use petgraph::graph::EdgeReference;
use petgraph::graph::{Graph, UnGraph, EdgeIndex};
use petgraph::graph::NodeIndex;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::Write;
use petgraph::visit::EdgeRef;
use std::collections::{HashSet, BTreeMap};

mod parse;
pub use parse::*;

mod visualize;
pub use visualize::*;

mod naming;
pub use naming::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    C,
    CAromatic,
    H,
    O,
    OAromatic,
    N,
    NAromatic,
    F,
    Cl,
    Br,
    // ... other elements
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Bond {
    Single,
    Double,
    Triple,
    Aromatic,
}

pub type MoleculeGraph = petgraph::graph::UnGraph<Element, Bond>;
