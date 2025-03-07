// use crate::*;
// use std::collections::HashMap;

// use petgraph::graph::{UnGraph, NodeIndex};
// use tracing::*;

// /// A fragment of a molecule
// #[derive(Debug, Clone)]
// pub struct Substituent {
//     attachment_point: NodeIndex,
//     attachment_bond: Bond,
//     component: MoleculeGraph
// }

// impl Substituent {
//     pub fn is_same_as(&self, other: &Substituent) -> bool {
//         // Perform a graph isomorphism check
        
//         // First, attach them to a parent molecule
//         let mut parent = MoleculeGraph::default();
//         parent.add_node(Element::H);
//         self.attach_to(&mut parent, NodeIndex::new(0));

//         let mut other_parent = MoleculeGraph::default();
//         parent.add_node(Element::H);
//         other.attach_to(&mut other_parent, NodeIndex::new(0));

//         // Then, check if the two graphs are isomorphic
//         petgraph::algo::is_isomorphic_matching(&parent, &other_parent, |a, b| a == b, |a, b| a == b)
//     }

//     pub fn root(&self) -> Element {
//         self.component[self.attachment_point]
//     }

//     pub fn amine() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::N);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn carboxyl() -> [Self; 2] {
//         [
//             Substituent::carbonyl(),
//             Substituent::hydroxy()
//         ]
//     }

//     pub fn aldehyde() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         Self { attachment_point, attachment_bond: Bond::Double, component: mol }
//     }

//     /// An -OH group
//     pub fn hydroxy() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn carbonyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         Self { attachment_point, attachment_bond: Bond::Double, component: mol }
//     }

//     pub fn methoxy() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         let carbon = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, carbon, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn ethoxy() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         let carbon1 = mol.add_node(Element::C);
//         let carbon2 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, carbon1, Bond::Single);
//         mol.add_edge(carbon1, carbon2, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn propoxy() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::O);
//         let carbon1 = mol.add_node(Element::C);
//         let carbon2 = mol.add_node(Element::C);
//         let carbon3 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, carbon1, Bond::Single);
//         mol.add_edge(carbon1, carbon2, Bond::Single);
//         mol.add_edge(carbon1, carbon3, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn acetyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let carbon1 = mol.add_node(Element::C);
//         let oxygen = mol.add_node(Element::O);
//         mol.add_edge(attachment_point, carbon1, Bond::Single);
//         mol.add_edge(carbon1, oxygen, Bond::Double);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn benzyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let carbon1 = mol.add_node(Element::C);
//         let carbon2 = mol.add_node(Element::C);
//         let carbon3 = mol.add_node(Element::C);
//         let carbon4 = mol.add_node(Element::C);
//         let carbon5 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, carbon1, Bond::Single);
//         mol.add_edge(carbon1, carbon2, Bond::Single);
//         mol.add_edge(carbon2, carbon3, Bond::Single);
//         mol.add_edge(carbon3, carbon4, Bond::Single);
//         mol.add_edge(carbon4, carbon5, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn phenyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let carbon1 = mol.add_node(Element::C);
//         let carbon2 = mol.add_node(Element::C);
//         let carbon3 = mol.add_node(Element::C);
//         let carbon4 = mol.add_node(Element::C);
//         let carbon5 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, carbon1, Bond::Aromatic);
//         mol.add_edge(carbon1, carbon2, Bond::Aromatic);
//         mol.add_edge(carbon2, carbon3, Bond::Aromatic);
//         mol.add_edge(carbon3, carbon4, Bond::Aromatic);
//         mol.add_edge(carbon4, carbon5, Bond::Aromatic);
//         mol.add_edge(carbon5, attachment_point, Bond::Aromatic);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn halide(halide: Element) -> Self {
//         assert!(halide.is_halogen());
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(halide);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn propyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);

//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(c1, c2, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn alkane(n: usize) -> Self {
//         Self::chain(vec![(Bond::Single, Element::C); n])
//     }

//     pub fn cyclic_alkane(n: usize) -> Self {
//         Self::ring(vec![(Bond::Single, Element::C); n])
//     }
    
//     pub fn isopropyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(attachment_point, c2, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn butyl() -> Self {
//         Self::alkane(4)
//     }

//     pub fn isobutyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);
//         let c3 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(c1, c2, Bond::Single);
//         mol.add_edge(c1, c3, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn sec_butyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);
//         let c3 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(attachment_point, c2, Bond::Single);
//         mol.add_edge(c2, c3, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn tert_butyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);
//         let c3 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(attachment_point, c2, Bond::Single);
//         mol.add_edge(attachment_point, c3, Bond::Single);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn vinyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Double);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn allyl() -> Self {
//         let mut mol = MoleculeGraph::default();
//         let attachment_point = mol.add_node(Element::C);
//         let c1 = mol.add_node(Element::C);
//         let c2 = mol.add_node(Element::C);
//         mol.add_edge(attachment_point, c1, Bond::Single);
//         mol.add_edge(c1, c2, Bond::Double);
//         Self { attachment_point, attachment_bond: Bond::Single, component: mol }
//     }

//     pub fn ring(atoms: Vec<(Bond, Element)>) -> Self {
//         let mut mol = MoleculeGraph::default();
//         let mut nodes = Vec::new();
//         let attachment_bond = atoms[0].0;
//         for (_, element) in &atoms {
//             let node = mol.add_node(*element);
//             nodes.push(node);
//         }
        

//         for i in 0..nodes.len() {
//             let j = (i + 1) % nodes.len();
//             mol.add_edge(nodes[i], nodes[j], atoms[i].0);
//         }

//         let attachment_point = nodes[0];

//         Self { attachment_point, attachment_bond, component: mol }
//     }

//     pub fn chain(atoms: Vec<(Bond, Element)>) -> Self {
//         let mut mol = MoleculeGraph::default();
//         let mut nodes = Vec::new();
//         let attachment_bond = atoms[0].0;
//         for (_, element) in &atoms {
//             let node = mol.add_node(*element);
//             nodes.push(node);
//         }
        

//         for i in 0..nodes.len() - 1 {
//             mol.add_edge(nodes[i], nodes[i + 1], atoms[i].0);
//         }

//         let attachment_point = nodes[0];

//         Self { attachment_point, attachment_bond, component: mol }
//     }

//     /// Fuse this substituent with another substituent at a given edge.
//     fn fuse(&self, child: Substituent, edge_idx_in_self: EdgeIndex, edge_idx_in_child: EdgeIndex) -> Self {
//         let mut mol = self.component.clone();
//         let other = child.component.clone();

//         let (self_source, self_target) = self.component.edge_endpoints(edge_idx_in_self).unwrap();
//         let (child_source, child_target) = child.component.edge_endpoints(edge_idx_in_child).unwrap();

//         // Assert that the sources are the same atom in both cases
//         assert_eq!(mol[self_source], other[child_source]);
//         // Assert that the targets are the same atom in both cases
//         assert_eq!(mol[self_target], other[child_target]);


//         // Add all the nodes from the child to the parent
//         let mut node_map = HashMap::new();
//         for node in other.node_indices() {
//             // Check if the node is a part of the edge
//             if node == child_source {
//                 node_map.insert(node, self_source);
//                 continue;
//             }

//             if node == child_target {
//                 node_map.insert(node, self_target);
//                 continue;
//             }

//             let new_node = mol.add_node(other[node]);
//             node_map.insert(node, new_node);
//         }

//         // Add all the edges from the child to the parent
//         for edge in other.edge_indices() {
//             if edge == edge_idx_in_child {
//                 continue;
//             }

//             let (source, target) = other.edge_endpoints(edge).unwrap();
//             trace!("Edge: {:?} -> {:?}", source, target);
//             let source = node_map[&source];
//             let target = node_map[&target];
//             trace!("Mapped edge: {:?} -> {:?}", source, target);
//             let bond = other[edge];
//             mol.add_edge(source, target, bond);
//         }

//         Self { attachment_point: self.attachment_point, attachment_bond: self.attachment_bond, component: mol }
//     }

//     fn nth_bond(&self, n: usize) -> Option<EdgeIndex> {
//         let mut count = 0;
//         for edge in self.component.edge_indices() {
//             if count == n {
//                 return Some(edge);
//             }
//             count += 1;
//         }
//         None
//     }

//     /// Attach this substituent to a parent molecule at a given node in the parent molecule.
//     fn attach_to(&self, parent: &mut MoleculeGraph, at: NodeIndex) {
//         let mut node_map = HashMap::new();
//         let attachment_point = parent.add_node(self.component[self.attachment_point]);
//         node_map.insert(self.attachment_point, attachment_point);

//         for node in self.component.node_indices() {
//             if node == self.attachment_point {
//                 continue;
//             }

//             let new_node = parent.add_node(self.component[node]);
//             node_map.insert(node, new_node);
//         }

//         for edge in self.component.edge_indices() {
//             let (source, target) = self.component.edge_endpoints(edge).unwrap();
//             let source = node_map[&source];
//             let target = node_map[&target];
//             let bond = self.component[edge];
//             parent.add_edge(source, target, bond);
//         }

//         parent.add_edge(at, attachment_point, self.attachment_bond);
//     }

//     pub fn fuse_ring(&self, other_ring: Substituent, nth_self: usize, nth_other: usize) -> Self {
//         let edge_idx_in_self = self.nth_bond(nth_self).unwrap();
//         let edge_idx_in_other = other_ring.nth_bond(nth_other).unwrap();
//         self.fuse(other_ring, edge_idx_in_self, edge_idx_in_other)
//     }

//     pub fn visualize(&self, filename: &str) {
//         visualize_graph(&self.component, "tmp.dot", Some(filename)).unwrap();
//     }
// }

// #[derive(Debug, Clone)]
// pub enum ParentChain {
//     Chain(Vec<(Bond, Element)>),
//     Ring(Vec<(Bond, Element)>),
// }

// impl ParentChain {
//     pub fn chain(atoms: Vec<(Bond, Element)>) -> Self {
//         Self::Chain(atoms)
//     }

//     pub fn ring(atoms: Vec<(Bond, Element)>) -> Self {
//         Self::Ring(atoms)
//     }

//     pub fn is_ring(&self) -> bool {
//         matches!(self, Self::Ring(_))
//     }

//     pub fn alkane(n: usize) -> Self {
//         Self::chain(vec![(Bond::Single, Element::C); n])
//     }

//     pub fn cyclic_alkane(n: usize) -> Self {
//         Self::ring(vec![(Bond::Single, Element::C); n])
//     }

//     fn iter(&self) -> impl Iterator<Item = &(Bond, Element)> {
//         match self {
//             Self::Chain(atoms) => atoms.iter(),
//             Self::Ring(atoms) => atoms.iter(),
//         }
//     }

//     fn to_graph(&self) -> (MoleculeGraph, Vec<NodeIndex>) {
//         let mut graph = MoleculeGraph::default();
//         let mut nodes = Vec::new();

//         for (i, (bond, element)) in self.iter().enumerate() {
//             let node = graph.add_node(*element);
//             nodes.push(node);

//             if i > 0 {
//                 graph.add_edge(nodes[i - 1], node, *bond);
//             }
//         }

//         if self.is_ring() {
//             let first = nodes[0];
//             let last = nodes[nodes.len() - 1];
//             graph.add_edge(first, last, self.iter().next().unwrap().0);
//         }

//         (graph, nodes)
//     }

//     /// Add the chain to the parent graph, attached at the `at` node.
//     /// 
//     /// Return a list of node indices in the parent graph that represent the Nth
//     /// atom in this parent chain.
//     fn attach_to(&self, parent: &mut MoleculeGraph, at: NodeIndex) -> Vec<NodeIndex> {
//         let mut nodes = Vec::new();

//         for (i, (bond, element)) in self.iter().enumerate() {
//             let node = parent.add_node(*element);
//             nodes.push(node);

//             if i == 0 {
//                 parent.add_edge(at, node, *bond);
//             } else {
//                 parent.add_edge(nodes[i - 1], node, *bond);
//             }
//         }

//         if self.is_ring() {
//             let first = nodes[0];
//             let last = nodes[nodes.len() - 1];
//             parent.add_edge(first, last, self.iter().next().unwrap().0);
//         }

//         nodes
//     }
// }

// impl From<Element> for ParentChain {
//     fn from(element: Element) -> Self {
//         ParentChain::chain(vec![(Bond::Single, element)])
//     }
// }

// #[derive(Debug, Clone)]
// pub struct FunctionalGroup {
//     parent_chain: ParentChain,
//     substituents: Vec<(usize, Substituent)>,
// }

// impl FunctionalGroup {
//     pub fn new(parent_chain: ParentChain) -> Self {
//         Self { parent_chain: parent_chain, substituents: Vec::new() }
//     }

//     pub fn alkane(n: usize) -> Self {
//         Self::new(ParentChain::alkane(n))
//     }

//     pub fn cyclic_alkane(n: usize) -> Self {
//         Self::new(ParentChain::cyclic_alkane(n))
//     }

//     pub fn attach_substituent(&mut self, at: usize, substituent: Substituent) {
//         self.substituents.push((at, substituent));
//     }

//     pub fn with_substituent(mut self, at: usize, substituent: Substituent) -> Self {
//         self.attach_substituent(at, substituent);
//         self
//     }

//     pub fn with_substituents(mut self, at: usize, substituents: impl IntoIterator<Item=Substituent>) -> Self {
//         for substituent in substituents {
//             self.attach_substituent(at, substituent);
//         }
//         self
//     }

//     pub fn to_graph(&self) -> MoleculeGraph {
//         let (mut graph, node_indices)  = self.parent_chain.to_graph();

//         trace!("Parent chain: {:#?}", graph);
//         for (at, substituent) in &self.substituents {
//             substituent.attach_to(&mut graph, node_indices[*at]);
//         }

//         graph
//     }

//     fn attach_to(&self, parent: &mut MoleculeGraph, at: NodeIndex) {
//         let node_indices = self.parent_chain.attach_to(parent, at);
//         for (at, substituent) in &self.substituents {
//             substituent.attach_to(parent, node_indices[*at]);
//         }
//     }

//     pub fn is_same_as(&self, other: &Self) -> bool {
//         // Construct the graph for the other molecule
//         let other_graph = other.to_graph();
//         let this_graph = self.to_graph();

//         // Check if the two graphs are isomorphic
//         petgraph::algo::is_isomorphic_matching(&this_graph, &other_graph, |a, b| a == b, |a, b| a == b)
//     }

//     pub fn visualize(&self, filename: &str) {
//         let mut graph = self.to_graph();
//         graph.hydrogenate();
//         visualize_graph(&graph, "tmp.dot", Some(filename)).unwrap();
//     }
// }

// pub enum OrganylGroup {
//     Ester {
//         acyl: Box<OrganylGroup>,
//         alkyl: Box<OrganylGroup>,
//     },
//     Ether {
//         alkyl1: Box<OrganylGroup>,
//         alkyl2: Box<OrganylGroup>,
//     },
//     Basic(FunctionalGroup),
// }

// impl OrganylGroup {
//     pub fn to_graph(&self) -> MoleculeGraph {
//         match self {
//             Self::Ester { acyl, alkyl } => {
//                 // let mut graph = acyl.to_graph();
//                 // alkyl.attach_to(&mut graph, NodeIndex::new(0));
//                 // graph

//                 let mut carboxyl = MoleculeGraph::default();
//                 let carbon = carboxyl.add_node(Element::C);
//                 let oxygen1 = carboxyl.add_node(Element::O);
//                 carboxyl.add_edge(carbon, oxygen1, Bond::Double);
//                 let oxygen2 = carboxyl.add_node(Element::O);
//                 carboxyl.add_edge(carbon, oxygen2, Bond::Single);

//                 acyl.attach_to(&mut carboxyl, oxygen1);
//                 alkyl.attach_to(&mut carboxyl, oxygen2);
//                 carboxyl

//             },
//             Self::Ether { alkyl1, alkyl2 } => {
//                 // let mut graph = alkyl1.to_graph();
//                 // alkyl2.attach_to(&mut graph, NodeIndex::new(0));
//                 // graph
//             },
//             Self::Basic(group) => group.to_graph(),
//         }
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_ring() {
//         let ring = Substituent::ring(vec![(Bond::Aromatic, Element::C); 6]);
//         ring.visualize("ring.png");
//     }

//     #[test]
//     fn test_fused_ring() {
//         init_logging(tracing::metadata::LevelFilter::TRACE);
//         let ring1 = Substituent::ring(vec![
//             Element::N, Element::C, Element::C, Element::C, Element::C, Element::C
//         ].iter().map(|&e| (Bond::Single, e)).collect());

//         let ring2 = Substituent::ring(vec![
//             Element::O, Element::C, Element::C, Element::C, Element::C, Element::C
//         ].iter().map(|&e| (Bond::Single, e)).collect());

//         let fused = ring1.fuse_ring(ring2, 1, 1);
//         fused.visualize("fused_ring.png");
//     }

//     #[test]
//     fn test_named_substituents() {
//         let base = FunctionalGroup::new(ParentChain::chain(vec![(Bond::Single, Element::O), (Bond::Single, Element::N)]));

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::propyl());
//         mol.visualize("ON-propyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::isopropyl());
//         mol.visualize("ON-isopropyl.png");


//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::butyl());
//         mol.visualize("ON-butyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::isobutyl());
//         mol.visualize("ON-isobutyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::sec_butyl());
//         mol.visualize("ON-sec-butyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::tert_butyl());
//         mol.visualize("ON-tert-butyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::vinyl());
//         mol.visualize("ON-vinyl.png");

//         let mut mol = base.clone();
//         mol.attach_substituent(1, Substituent::allyl());
//         mol.visualize("ON-allyl.png");
//     }

//     #[test]
//     fn test_alkanes() {
//         init_logging(tracing::metadata::LevelFilter::TRACE);

//         let hexane = FunctionalGroup::alkane(6);
//         hexane.visualize("hexane.png");

//         assert!(!hexane.is_same_as(&FunctionalGroup::alkane(5)));
//         assert!(hexane.is_same_as(&FunctionalGroup::alkane(6)));
//         assert!(!hexane.is_same_as(&FunctionalGroup::cyclic_alkane(6)));
//         assert!(!hexane.is_same_as(&FunctionalGroup::alkane(7)));
//     }

//     #[test]
//     fn test_named_substituents2() {
//         let ethanol = FunctionalGroup::alkane(2)
//             .with_substituent(1, Substituent::hydroxy());
//         ethanol.visualize("ethanol.png");

//         let phenylethanol = FunctionalGroup::alkane(2)
//             .with_substituent(0, Substituent::phenyl())
//             .with_substituent(1, Substituent::hydroxy());
//         phenylethanol.visualize("phenylethanol.png");

//         let acetone = FunctionalGroup::alkane(3)
//             .with_substituent(1, Substituent::carbonyl());
//         acetone.visualize("acetone.png");

//         let phenylalanine = FunctionalGroup::alkane(3)
//             .with_substituent(0, Substituent::phenyl())
//             .with_substituent(1, Substituent::amine())
//             .with_substituents(2, Substituent::carboxyl());
//         phenylalanine.visualize("phenylalanine.png");

//         let ethanoic_acid = FunctionalGroup::alkane(2)
//             .with_substituents(1, Substituent::carboxyl());
//         ethanoic_acid.visualize("acetic-acid.png");

//         let formaldehyde = FunctionalGroup::alkane(1)
//             .with_substituent(0, Substituent::aldehyde());
//         formaldehyde.visualize("formaldehyde.png");
//     }
// }