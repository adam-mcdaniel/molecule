
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

pub fn init_logging(level: tracing::metadata::LevelFilter) {
    let _ = tracing_subscriber::fmt()
        .with_max_level(level)
        .without_time()
        // .with_target(false)
        .try_init();
}

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
    S,
    SAromatic
    // ... other elements
}

impl Element {
    pub fn strip_aromatic(&self) -> Element {
        match self {
            Element::CAromatic => Element::C,
            Element::OAromatic => Element::O,
            Element::NAromatic => Element::N,
            Element::SAromatic => Element::S,
            _ => *self,
        }
    }

    pub fn is_carbon(&self) -> bool {
        match self {
            Element::C | Element::CAromatic => true,
            _ => false,
        }
    }

    pub fn is_hydrogen(&self) -> bool {
        match self {
            Element::H => true,
            _ => false,
        }
    }

    pub fn is_oxygen(&self) -> bool {
        match self {
            Element::O | Element::OAromatic => true,
            _ => false,
        }
    }

    pub fn is_nitrogen(&self) -> bool {
        match self {
            Element::N | Element::NAromatic => true,
            _ => false,
        }
    }

    pub fn is_sulfur(&self) -> bool {
        match self {
            Element::S | Element::SAromatic => true,
            _ => false,
        }
    }

    pub fn symbol(&self) -> &'static str {
        match self {
            Element::C => "C",
            Element::CAromatic => "C",
            Element::H => "H",
            Element::O => "O",
            Element::OAromatic => "O",
            Element::N => "N",
            Element::NAromatic => "N",
            Element::F => "F",
            Element::Cl => "Cl",
            Element::Br => "Br",
            Element::S => "S",
            Element::SAromatic => "S",
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

pub type MoleculeGraph = petgraph::graph::UnGraph<Element, Bond>;
