// src/parser.rs

use crate::{Bond, Bond::*, Element, Element::*};
use crate::*;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use crate::MoleculeGraph;

// src/parser.rs
use std::fmt::Write; // Correctly import the Write trait

/// Helper function to determine if an element is aromatic
fn is_aromatic(element: &Element) -> bool {
    matches!(element, CAromatic | NAromatic | OAromatic | SAromatic)
}

/// Parses a SMILES string into a MoleculeGraph.
///
/// # Arguments
///
/// * `smiles` - The SMILES string to parse.
///
/// # Returns
///
/// * `Result<MoleculeGraph, String>` - The parsed molecular graph or an error message.
pub fn parse_smiles(smiles: &str) -> Result<MoleculeGraph, String> {
    let mut graph = MoleculeGraph::new_undirected();
    let mut current_atom: Option<NodeIndex> = None;
    let mut bond_type = Single; // Default bond type
    let mut branch_stack: Vec<NodeIndex> = Vec::new();
    let mut ring_map: HashMap<u8, NodeIndex> = HashMap::new();

    let chars: Vec<char> = smiles.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        let c = chars[i];

        match c {
            '(' => {
                // Start of a branch: push current_atom to stack
                if let Some(atom) = current_atom {
                    branch_stack.push(atom);
                } else {
                    return Err(format!(
                        "Branch start '(' at position {} without a current atom",
                        i
                    ));
                }
                i += 1;
            }
            ')' => {
                // End of a branch: pop from stack
                current_atom = branch_stack.pop();
                if current_atom.is_none() {
                    return Err(format!(
                        "Branch end ')' at position {} without a matching '('",
                        i
                    ));
                }
                i += 1;
            }
            '-' | '=' | '#' | ':' => {
                // Bond symbol
                bond_type = match c {
                    '-' => Single,
                    '=' => Double,
                    '#' => Triple,
                    ':' => Aromatic, // Assume ':' represents aromatic bond
                    _ => {
                        return Err(format!(
                            "Unknown bond type '{}' at position {}",
                            c, i
                        ))
                    }
                };
                i += 1;
            }
            '0'..='9' => {
                // Ring closure
                let ring_number = c.to_digit(10).unwrap() as u8;
                if let Some(&start_atom) = ring_map.get(&ring_number) {
                    // Create bond between current_atom and start_atom
                    let bond_to_use = if let Some(current) = current_atom {
                        if is_aromatic(&graph[start_atom]) && is_aromatic(&graph[current]) {
                            Bond::Aromatic
                        } else {
                            bond_type.clone()
                        }
                    } else {
                        bond_type.clone()
                    };
                    graph.add_edge(
                        current_atom.ok_or(format!(
                            "Ring closure digit '{}' at position {} without a current atom",
                            ring_number, i
                        ))?,
                        start_atom,
                        bond_to_use,
                    );
                    // Reset bond_type to default after use
                    bond_type = Single;
                    // Remove ring closure
                    ring_map.remove(&ring_number);
                } else {
                    // Store the current_atom for future ring closure
                    ring_map.insert(
                        ring_number,
                        current_atom.ok_or(format!(
                            "Ring closure digit '{}' at position {} without a current atom",
                            ring_number, i
                        ))?,
                    );
                }
                i += 1;
            }
            '[' => {
                // Bracketed atom (handling atom properties is beyond this basic parser)
                // For simplicity, we'll skip bracketed atoms in this example
                // You can extend this section to handle complex atoms
                // Find the closing bracket
                if let Some(end_relative) = chars[i..].iter().position(|&x| x == ']') {
                    // Calculate the absolute position of the closing bracket
                    let end = i + end_relative;
                    // Extract the content inside brackets
                    let atom_str: String = chars[i + 1..end].iter().collect();
                    let element = parse_element(&atom_str).expect("Failed to parse element in brackets");
                    let new_atom = graph.add_node(element);
                    if let Some(prev_atom) = current_atom {
                        let bond_to_use = if is_aromatic(&graph[prev_atom]) && is_aromatic(&element) {
                            Bond::Aromatic
                        } else {
                            bond_type.clone()
                        };
                        graph.add_edge(prev_atom, new_atom, bond_to_use);
                        bond_type = Single;
                    }
                    current_atom = Some(new_atom);
                    i = end + 1;
                } else {
                    return Err(format!("Unclosed bracket '[' at position {}", i));
                }
            }
            _ => {

                let mut atom_str = String::new();
                // Check if there's a next character that is lowercase.
                if i + 1 < chars.len() && chars[i+1].is_lowercase() {
                    // Form the two-character candidate.
                    let candidate = &smiles[i..i+2];
                    // If the candidate forms a valid element...
                    if parse_element(candidate).is_ok() {
                        // If current char is uppercase, always use the two-letter element.
                        if c.is_uppercase() {
                            atom_str.push_str(candidate);
                            i += 2;
                        } else {
                            // For aromatic (lowercase) atoms, check if the second char by itself is valid.
                            let second = &smiles[i+1..i+2];
                            if parse_element(second).is_ok() {
                                // The second letter is valid on its own, so treat current as a single atom.
                                atom_str.push(c);
                                i += 1;
                            } else {
                                // Otherwise, combine them (e.g. "cl" for aromatic chlorine).
                                atom_str.push_str(candidate);
                                i += 2;
                            }
                        }
                    } else {
                        // If the two-letter candidate isn’t valid, just take the single letter.
                        atom_str.push(c);
                        i += 1;
                    }
                } else {
                   // No next lowercase char available; use the current character alone.
                   atom_str.push(c);
                   i += 1;
                }
            
                // Now parse the element from atom_str.
                let element = parse_element(&atom_str).expect("Failed to parse element");
                let new_atom = graph.add_node(element);
                if let Some(prev_atom) = current_atom {
                    let is_prev_aromatic = is_aromatic(&graph[prev_atom]);
                    let is_new_aromatic = is_aromatic(&element);
                    let bond_to_use = if is_prev_aromatic && is_new_aromatic {
                        Bond::Aromatic
                    } else {
                        bond_type.clone()
                    };
                    graph.add_edge(prev_atom, new_atom, bond_to_use);
                    bond_type = Single; // Reset bond type after use.
                }
                current_atom = Some(new_atom);
            }
        }
    }

    Ok(graph)
}

fn parse_element(symbol: &str) -> Result<Element, String> {
    match symbol {
        "C" => Ok(Element::C),
        "c" => Ok(Element::CAromatic),
        "H" => Ok(Element::H),
        "O" => Ok(Element::O),
        "o" => Ok(Element::OAromatic),
        "N" => Ok(Element::N),
        "n" => Ok(Element::NAromatic),
        "Cl" => Ok(Element::Cl),
        "Br" => Ok(Element::Br),
        "F" => Ok(Element::F),
        "S" => Ok(Element::S),
        "s" => Ok(Element::SAromatic),
        // ... handle other elements as needed
        _ => Err(format!("Unknown element symbol '{}'", symbol)),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;

    // Assuming visualize_graph is implemented as per the previous assistant's message
    // Ensure you have the visualize_graph function in your parser.rs

    #[test]
    fn test_parse_ethanol() {
        let smiles = "CCO"; // Ethanol
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "ethanol.dot", Some("ethanol.png"))
            .expect("Failed to visualize graph");

        assert_eq!(molecule.node_count(), 3); // 2 Carbons, 1 Oxygen

        // Check elements
        assert_eq!(molecule[NodeIndex::new(0)], C);
        assert_eq!(molecule[NodeIndex::new(1)], C);
        assert_eq!(molecule[NodeIndex::new(2)], O);

        // Check bonds
        let edges: Vec<_> = molecule.edge_references().collect();
        assert_eq!(edges.len(), 2); // C-C and C-O

        for edge in edges {
            let source = edge.source().index();
            let target = edge.target().index();
            let bond = edge.weight();

            match (source, target) {
                (0, 1) | (1, 0) => assert_eq!(bond, &Single),
                (1, 2) | (2, 1) => assert_eq!(bond, &Single),
                _ => panic!("Unexpected bond between {:?} and {:?}", source, target),
            }
        }
    }

    #[test]
    fn test_parse_cyclohexane() {
        let smiles = "C1CCCCC1"; // Cyclohexane
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "cyclohexane.dot", Some("cyclohexane.png"))
            .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
        // assert_eq!(molecule.node_count(), 6); // 6 Carbons

        // Check all atoms are carbon
        for node in molecule.node_indices() {
            assert_eq!(molecule[node], C);
        }

        // Check bonds
        let edges: Vec<_> = molecule.edge_references().collect();
        assert_eq!(edges.len(), 6); // C-C bonds forming a ring

        // Each carbon should have two bonds (since it's a ring)
        for node in molecule.node_indices() {
            let degree = molecule.edges(node).count();
            assert_eq!(degree, 2, "Node {} has degree {}", node.index(), degree);
        }
        println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_ciprofloxacin() {
        let smiles = "C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O"
            .to_string(); // Ciprofloxacin

        let molecule = parse_smiles(&smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "ciprofloxacin.dot", Some("ciprofloxacin.png"))
            .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
        println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_1_choloro_ethane() {
        let smiles = "CCCl"
            .to_string(); // Ciprofloxacin

        let molecule = parse_smiles(&smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "1-choloro-ethane.dot", Some("1-choloro-ethane.png"))
            .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
        println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_isobutane() {
        let smiles = "CC(C)C"; // Isobutane
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "isobutane.dot", Some("isobutane.png"))
            .expect("Failed to visualize graph");

        assert_eq!(molecule.node_count(), 4); // 4 Carbons

        // Check elements
        for node in molecule.node_indices() {
            assert_eq!(molecule[node], C);
        }

        // Check bonds
        let edges: Vec<_> = molecule.edge_references().collect();
        assert_eq!(edges.len(), 3); // C-C bonds

        // In an undirected graph, cyclohexane has 6 edges, isobutane should have 3 unique edges
        assert_eq!(edges.len(), 3, "Isobutane should have 3 unique C-C bonds");

        // Verify specific bonds
        let bonds: Vec<_> = molecule.edge_references().collect();

        // Expected bonds:
        // C0 - C1
        // C1 - C2
        // C1 - C3

        assert!(bonds.iter().any(|e| {
            (e.source() == NodeIndex::new(0) && e.target() == NodeIndex::new(1))
                || (e.source() == NodeIndex::new(1) && e.target() == NodeIndex::new(0))
        }));

        assert!(bonds.iter().any(|e| {
            (e.source() == NodeIndex::new(1) && e.target() == NodeIndex::new(2))
                || (e.source() == NodeIndex::new(2) && e.target() == NodeIndex::new(1))
        }));

        assert!(bonds.iter().any(|e| {
            (e.source() == NodeIndex::new(1) && e.target() == NodeIndex::new(3))
                || (e.source() == NodeIndex::new(3) && e.target() == NodeIndex::new(1))
        }));
        println!("Name: {}", iupac_name(&molecule));
    }
}