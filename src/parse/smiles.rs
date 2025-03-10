// src/parser.rs

use crate::MoleculeGraph;
use crate::*;
use crate::{Bond, Bond::*, Element};
use anyhow::Result;
use petgraph::graph::NodeIndex;
use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

// src/parser.rs
use std::fmt::Write; // Correctly import the Write trait

#[derive(Error, Debug)]
pub enum SmilesError {
    #[error("Branch start '(' at position {0} (followed by {1}) without a current atom")]
    BranchNoCurrentAtom(usize, String),
    #[error("Branch end ')' at position {0} (followed by {1}) without a matching '('")]
    BranchEndNoStart(usize, String),
    #[error("Ring closure digit '{0}' at position {1} without a current atom")]
    RingClosureNoCurrentAtom(char, usize),
    #[error("Unclosed bracket '[' at position {0}")]
    UnclosedBracket(usize),
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
pub fn parse_smiles(smiles: &str) -> Result<MoleculeGraph> {
    parse_smiles_helper(smiles)
        .context(format!("Failed to parse SMILES string {smiles}"))
}

/// Helper: Merge a fragment MoleculeGraph into the main graph.
/// Returns a mapping from the fragment’s NodeIndex to the new NodeIndex in `main`.
fn merge_fragment_into_graph(main: &mut MoleculeGraph, frag: &MoleculeGraph) -> BTreeMap<NodeIndex, NodeIndex> {
    let mut mapping = BTreeMap::new();
    // Add all nodes from the fragment into the main graph.
    for node in frag.node_indices() {
        // Clone the element (assuming Element implements Clone).
        let element = frag[node].clone();
        let new_node = main.add_node(element);
        mapping.insert(node, new_node);
    }
    // Add all edges from the fragment into the main graph.
    for edge in frag.edge_references() {
        let src = mapping[&edge.source()];
        let tgt = mapping[&edge.target()];
        main.add_edge(src, tgt, edge.weight().clone());
    }
    mapping
}

fn parse_bracketed_smiles(content: &str) -> anyhow::Result<(Element, usize)> {
    // Example content: "PH3", "Cl", "C", etc.
    // First, extract the element symbol.
    let mut chars = content.chars().peekable();
    let mut symbol = String::new();
    if let Some(c) = chars.next() {
        if c.is_uppercase() {
            symbol.push(c);
            // Try to take a second lowercase letter if available.
            if let Some(&next) = chars.peek() {
                if next.is_lowercase() {
                    symbol.push(chars.next().unwrap());
                }
            }
        } else {
            return Err(anyhow::anyhow!(
                "Bracket content must start with an uppercase letter, got {}",
                content
            ));
        }
    } else {
        return Err(anyhow::anyhow!("Empty bracket content"));
    }
    // Now, check if the next character is 'H'. If so, parse an explicit hydrogen count.
    let mut h_count = 0;
    if let Some(&next) = chars.peek() {
        if next == 'H' {
            // Consume the 'H'
            chars.next();
            let mut digit_str = String::new();
            while let Some(&d) = chars.peek() {
                if d.is_digit(10) {
                    digit_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }
            h_count = if digit_str.is_empty() {
                1
            } else {
                digit_str.parse::<usize>()?
            };
        }
    }
    let element = Element::from_smiles(&symbol)
        .context(format!("Failed to parse element symbol: {}", symbol))?;
    Ok((element, h_count))
}

fn parse_smiles_helper(smiles: &str) -> Result<MoleculeGraph> {
    let mut graph = MoleculeGraph::new_undirected();
    let mut current_atom: Option<NodeIndex> = None;
    let mut bond_type = Single; // Default bond type
    let mut branch_stack: Vec<NodeIndex> = Vec::new();
    let mut ring_map: BTreeMap<u8, NodeIndex> = BTreeMap::new();

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
                    return Err(SmilesError::BranchNoCurrentAtom(i, smiles[i..].to_string()))
                        .context(format!("While parsing {smiles}"));
                }
                i += 1;
            }
            ')' => {
                // End of a branch: pop from stack
                current_atom = branch_stack.pop();
                if current_atom.is_none() {
                    // return Err(format!(
                    //     "Branch end ')' at position {} without a matching '('",
                    //     i
                    // ));
                    return Err(SmilesError::BranchEndNoStart(i, smiles[i..].to_string()))
                        .context(format!("While parsing {smiles}"));
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
                        return Err(SmilesError::BranchEndNoStart(i, smiles[i..].to_string()))
                            .context(format!("While parsing {smiles}"));
                    }
                };
                i += 1;
            }
            '%' => {
                // Ring closure with two-digit label.
                // Ensure there are two more characters.
                if i + 2 >= chars.len() {
                    return Err(anyhow::anyhow!("Incomplete ring closure after '%' at position {}", i));
                }
                let digits: String = chars[i+1..i+3].iter().collect();
                let ring_number: u8 = digits.parse().map_err(|e| {
                    anyhow::anyhow!("Failed to parse ring closure digits '{}' at position {}: {}", digits, i, e)
                })?;
                // Process the ring closure similarly to a single-digit closure.
                if let Some(&start_atom) = ring_map.get(&ring_number) {
                    let bond_to_use = if let Some(current) = current_atom {
                        if graph[start_atom].is_aromatic() && graph[current].is_aromatic() {
                            Bond::Aromatic
                        } else {
                            bond_type.clone()
                        }
                    } else {
                        bond_type.clone()
                    };
                    graph.add_edge(
                        current_atom
                            .ok_or(SmilesError::RingClosureNoCurrentAtom(c, i))
                            .context(format!("While parsing {smiles}"))?,
                        start_atom,
                        bond_to_use,
                    );
                    bond_type = Single;
                    ring_map.remove(&ring_number);
                } else {
                    ring_map.insert(
                        ring_number,
                        current_atom
                            .ok_or(SmilesError::RingClosureNoCurrentAtom(c, i))
                            .context(format!("While parsing {smiles}"))?,
                    );
                }
                i += 3; // Skip '%' and the two digits
            }
            '0'..='9' => {
                // Ring closure
                let ring_number = c.to_digit(10).unwrap() as u8;
                if let Some(&start_atom) = ring_map.get(&ring_number) {
                    // Create bond between current_atom and start_atom
                    let bond_to_use = if let Some(current) = current_atom {
                        if graph[start_atom].is_aromatic() && graph[current].is_aromatic() {
                            Bond::Aromatic
                        } else {
                            bond_type.clone()
                        }
                    } else {
                        bond_type.clone()
                    };
                    graph.add_edge(
                        // current_atom.ok_or(format!(
                        //     "Ring closure digit '{}' at position {} without a current atom",
                        //     ring_number, i
                        // ))?,
                        current_atom
                            .ok_or(SmilesError::RingClosureNoCurrentAtom(c, i))
                            .context(format!("While parsing {smiles}"))?,
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
                        // current_atom.ok_or(format!(
                        //     "Ring closure digit '{}' at position {} without a current atom",
                        //     ring_number, i
                        // ))?,
                        current_atom
                            .ok_or(SmilesError::RingClosureNoCurrentAtom(c, i))
                            .context(format!("While parsing {smiles}"))?,
                    );
                }
                i += 1;
            }
            '[' => {
                // Bracketed atom or molecule fragment.
                if let Some(end_relative) = chars[i..].iter().position(|&x| x == ']') {
                    let end = i + end_relative;
                    let raw_content: String = chars[i + 1..end].iter().collect();
                    // Remove any '+' or '-' characters (charges).
                    let cleaned: String = raw_content
                        .chars()
                        .filter(|&ch| ch != '+' && ch != '-')
                        .collect();
                    // Remove any leading digits (isotope information)
                    let cleaned = cleaned.trim_start_matches(|c: char| c.is_numeric()).to_string();

                    // Try to parse the cleaned content as an element with explicit hydrogen count.
                    if let Ok((element, h_count)) = parse_bracketed_smiles(&cleaned) {
                        let new_atom = graph.add_node(element);
                        if let Some(prev_atom) = current_atom {
                            let bond_to_use = if graph[prev_atom].is_aromatic() && graph[new_atom].is_aromatic() {
                                Bond::Aromatic
                            } else {
                                bond_type.clone()
                            };
                            graph.add_edge(prev_atom, new_atom, bond_to_use);
                            bond_type = Single;
                        }
                        current_atom = Some(new_atom);
                        // Add explicit hydrogens.
                        for _ in 0..h_count {
                            let h_node = graph.add_node(Element::H);
                            graph.add_edge(new_atom, h_node, Bond::Single);
                        }
                    } else {
                        // Otherwise, treat the cleaned content as a molecule fragment.
                        let fragment_graph = parse_smiles_helper(&cleaned)
                            .context(format!("Parsing fragment '{}' inside brackets", cleaned))?;
                        if let Some(prev_atom) = current_atom {
                            let mapping = merge_fragment_into_graph(&mut graph, &fragment_graph);
                            let frag_root = mapping.get(&NodeIndex::new(0))
                                .expect("Fragment has no nodes");
                            graph.add_edge(prev_atom, *frag_root, Single);
                            current_atom = Some(*frag_root);
                            bond_type = Single;
                        } else {
                            let mapping = merge_fragment_into_graph(&mut graph, &fragment_graph);
                            let frag_root = mapping.get(&NodeIndex::new(0))
                                .expect("Fragment has no nodes");
                            current_atom = Some(*frag_root);
                            bond_type = Single;
                        }
                    }
                    i = end + 1;
                } else {
                    return Err(SmilesError::UnclosedBracket(i))
                        .context(format!("While parsing {smiles}"));
                }
                /*
                // Bracketed atom (handling atom properties is beyond this basic parser)
                // For simplicity, we'll skip bracketed atoms in this example
                // You can extend this section to handle complex atoms
                // Find the closing bracket
                if let Some(end_relative) = chars[i..].iter().position(|&x| x == ']') {
                    // Calculate the absolute position of the closing bracket
                    let end = i + end_relative;
                    // Extract the content inside brackets
                    let atom_str: String = chars[i + 1..end].iter().collect();
                    let element = Element::from_smiles(&atom_str)
                        .context("Failed to parse element")?;
                    let new_atom = graph.add_node(element);
                    if let Some(prev_atom) = current_atom {
                        let bond_to_use = if graph[prev_atom].is_aromatic() && element.is_aromatic()
                        {
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
                    return Err(SmilesError::UnclosedBracket(i))
                        .context(format!("While parsing {smiles}"));
                }
                */
            }
            '@' | '/' | '\\' => {
                // Ignore stereochemistry markers, including slashes.
                i += 1;
            }
            '.' => {
                // Dot denotes that the next fragment is disconnected.
                // Reset the current atom (and branch stack) to start a new component.
                current_atom = None;
                branch_stack.clear();
                i += 1;
            }
            'R' => {
                // Count the number of apostrophes to determine the number of R groups.
                let mut r_count = 0;
                i += 1;
                while i < chars.len() && chars[i] == '\'' {
                    i += 1;
                    r_count += 1;
                }
                // Create a new atom for the R group.
                let atom_str = format!("R{}", "'".repeat(r_count));
                
                let element = Element::from_smiles(&atom_str).context("Failed to parse element")?;
                let new_atom = graph.add_node(element);
                if let Some(prev_atom) = current_atom {
                    let is_prev_aromatic = graph[prev_atom].is_aromatic();
                    let is_new_aromatic = element.is_aromatic();
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
            _ => {
                let mut atom_str = String::new();
                // Check if there's a next character that is lowercase.
                if i + 1 < chars.len() && chars[i + 1].is_lowercase() {
                    // Form the two-character candidate.
                    let candidate = &smiles[i..i + 2];
                    // If the candidate forms a valid element...
                    if Element::from_smiles(candidate).is_ok() {
                        // If current char is uppercase, always use the two-letter element.
                        if c.is_uppercase() {
                            atom_str.push_str(candidate);
                            i += 2;
                        } else {
                            // For aromatic (lowercase) atoms, check if the second char by itself is valid.
                            let second = &smiles[i + 1..i + 2];
                            if Element::from_smiles(second).is_ok() {
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
                let element = Element::from_smiles(&atom_str).context("Failed to parse element")?;
                let new_atom = graph.add_node(element);
                if let Some(prev_atom) = current_atom {
                    let is_prev_aromatic = graph[prev_atom].is_aromatic();
                    let is_new_aromatic = element.is_aromatic();
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


/// Convert a SMILES string into a canonical SMILES for database lookup.
/// 
/// This function uses your custom SMILES parser to get a molecule graph and then
/// performs a naïve DFS-based serialization (which works well only for acyclic graphs).
/// 
/// # Arguments
///
/// * `smiles` - A SMILES string to canonicalize.
///
/// # Returns
///
/// A `Result<String>` containing the canonical SMILES string or an error.
pub fn canonize_smiles(smiles: &str) -> Result<String> {
    init_logging("trace");
    // Parse the input into your MoleculeGraph.
    // let graph = parse_smiles(smiles)?;
    // // Compute a canonical SMILES from the graph.
    // let canon = compute_canonical_smiles(&graph);
    // Ok(canon)
    // Ok(smiles.to_string())
    let graph = parse_smiles(smiles)?;
    let graph = graph.morgan_canonize();
    // let graph = graph.canonize();
    let canon_smiles = molecule_to_smiles(&graph)
        .context(format!("Failed to canonize SMILES: {smiles}"))?;
    let original = Substituent::parse_smiles(smiles).expect("Failed to parse original SMILES");
    let canon = Substituent::parse_smiles(&canon_smiles).expect("Failed to parse canonical SMILES");
    if !original.is_same_as(&canon) {
        println!("Original: {}", smiles);
        println!("Canonical: {}", canon_smiles);

        original.visualize("original.png").expect("Failed to visualize original");
        canon.visualize("canonical.png").expect("Failed to visualize canonical");
        panic!("Original and canonical SMILES do not match");
    }
    Ok(canon_smiles)
}

/// A naïve DFS-based canonical SMILES generator.
/// 
/// **Note:** This implementation does not handle rings/cycles correctly. It is intended
/// only as an illustrative example of how you might “canonicalize” a molecule for lookup.
/// For a full solution, you would need to implement ring closure handling and a robust
/// canonical labeling algorithm.
fn compute_canonical_smiles(graph: &MoleculeGraph) -> String {
    // We choose a starting node. For simplicity, we use the node with index 0.
    let start = NodeIndex::new(0);
    let mut visited = BTreeSet::new();
    let mut result = String::new();
    dfs_smiles(graph, start, &mut visited, &mut result);
    result
}

/// A recursive DFS that appends node (atom) symbols and branches in a reproducible way.
///
/// For each node, we print its element symbol, then recursively process its unvisited
/// neighbors in sorted order (by the string representation of the neighbor’s element).
/// Branches are delimited by parentheses.
fn dfs_smiles(
    graph: &MoleculeGraph,
    node: NodeIndex,
    visited: &mut BTreeSet<NodeIndex>,
    result: &mut String,
) {
    visited.insert(node);
    // Write the element symbol for the current node.
    // (This assumes your Element type implements Display appropriately.)
    write!(result, "{}", graph[node].smiles_symbol()).expect("Write failed");

    // Collect all neighbor nodes that have not yet been visited.
    let mut neighbors: Vec<NodeIndex> = graph
        .neighbors(node)
        .filter(|n| !visited.contains(n))
        .collect();
    // Sort neighbors by their element symbol for reproducible order.
    neighbors.sort_by_key(|n| format!("{}", graph[*n].smiles_symbol()));

    // Process each neighbor as a branch.
    let mut remaining = neighbors.len();
    for neighbor in neighbors {
        if remaining == 1 {
            // If this is the last neighbor, we can append it directly.
            dfs_smiles(graph, neighbor, visited, result);
            continue;
        }
        result.push('(');
        dfs_smiles(graph, neighbor, visited, result);
        result.push(')');
        remaining -= 1;
    }
}

/// A helper structure to record a ring closure edge.
#[derive(Debug, Clone)]
struct RingClosure {
    opening: NodeIndex, // the ancestor where the ring is opened
    closing: NodeIndex, // the descendant where the ring is closed
    bond: Bond,         // the bond used for the closure
    digit: usize,       // assigned ring closure digit
}

// Helper function to check if a slice of bonds (of length 6) is alternating
fn is_alternating(bonds: &[Bond]) -> bool {
    if bonds.len() != 6 {
        return false;
    }
    let pattern1 = [
        Bond::Single, Bond::Double, Bond::Single, Bond::Double, Bond::Single, Bond::Double,
    ];
    let pattern2 = [
        Bond::Double, Bond::Single, Bond::Double, Bond::Single, Bond::Double, Bond::Single,
    ];
    bonds == pattern1 || bonds == pattern2
}

/// First pass: Compute a spanning tree of the graph (via DFS) and record ring closures.
/// Returns a mapping from each node (except the root) to its parent in the DFS tree,
/// and a vector of RingClosure records for every back edge encountered.
fn compute_spanning_tree_and_ring_closures(
    graph: &mut MoleculeGraph,
    start: NodeIndex,
) -> (BTreeMap<NodeIndex, NodeIndex>, Vec<RingClosure>) {
    let mut parent: BTreeMap<NodeIndex, NodeIndex> = BTreeMap::new();
    let mut ring_closures: Vec<RingClosure> = Vec::new();
    let mut visited: BTreeSet<NodeIndex> = BTreeSet::new();
    let mut path: Vec<NodeIndex> = Vec::new(); // current DFS stack

    fn dfs(
        graph: &mut MoleculeGraph,
        current: NodeIndex,
        parent: &mut BTreeMap<NodeIndex, NodeIndex>,
        ring_closures: &mut Vec<RingClosure>,
        visited: &mut BTreeSet<NodeIndex>,
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
                // Compute the cycle (ring) from nbr to current.
                if let Some(pos) = path.iter().position(|&node| node == nbr) {
                    let ring_nodes: Vec<NodeIndex> = path[pos..].to_vec(); // from ancestor to current
                    // The full cycle is ring_nodes plus the closure edge (current->nbr)
                    if ring_nodes.len() == 6 {
                        // Check that all atoms in the ring are carbon.
                        if ring_nodes.iter().all(|&n| {
                            // Assume Element has a method is_carbon()
                            graph[n].is_carbon()
                        }) {
                            // Collect the bonds along the ring.
                            let mut bonds = Vec::new();
                            for window in ring_nodes.windows(2) {
                                if let Some(edge) = graph.find_edge(window[0], window[1]) {
                                    bonds.push(*graph.edge_weight(edge).unwrap());
                                }
                            }
                            // Include the closure edge from current to nbr.
                            if let Some(edge) = graph.find_edge(*ring_nodes.last().unwrap(), nbr) {
                                bonds.push(*graph.edge_weight(edge).unwrap());
                            }
                            // Check if bonds alternate between single and double.
                            if is_alternating(&bonds) {
                                // Update all bonds in the cycle to aromatic.
                                for window in ring_nodes.windows(2) {
                                    if let Some(edge) = graph.find_edge(window[0], window[1]) {
                                        *graph.edge_weight_mut(edge).unwrap() = Bond::Aromatic;
                                    }
                                }
                                if let Some(edge) = graph.find_edge(*ring_nodes.last().unwrap(), nbr) {
                                    *graph.edge_weight_mut(edge).unwrap() = Bond::Aromatic;
                                }

                                // Update the elements in the ring to aromatic.
                                for &n in &ring_nodes {
                                    graph[n].set_aromatic(true);
                                }
                            }
                        }
                    }
                }
                // Record the ring closure if not already recorded.
                let key = if current < nbr { (current, nbr) } else { (nbr, current) };
                if !ring_closures.iter().any(|rc| {
                    let existing_key = if rc.opening < rc.closing {
                        (rc.opening, rc.closing)
                    } else {
                        (rc.closing, rc.opening)
                    };
                    existing_key == key
                }) {
                    // Record the closure with a placeholder digit.
                    if let Some(edge) = graph.find_edge(current, nbr) {
                        ring_closures.push(RingClosure {
                            opening: nbr,
                            closing: current,
                            bond: *graph.edge_weight(edge).unwrap(),
                            digit: 0,
                        });
                    }
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

// /// Second pass: Generate SMILES from the spanning tree, inserting ring–closure markers
// /// and forcing branch notation when the current atom is a branching center.
// fn generate_smiles_from_tree(
//     graph: &MoleculeGraph,
//     root: NodeIndex,
//     parent: &BTreeMap<NodeIndex, NodeIndex>,
//     ring_closures: &mut Vec<RingClosure>,
// ) -> String {
//     // First, assign digits to ring closures.
//     ring_closures.sort_by_key(|rc| (rc.opening.index(), rc.closing.index()));
//     for (i, rc) in ring_closures.iter_mut().enumerate() {
//         rc.digit = i + 1;
//     }

//     // Build the children mapping from the parent mapping.
//     let mut children: BTreeMap<NodeIndex, Vec<NodeIndex>> = BTreeMap::new();
//     for (&child, &p) in parent.iter() {
//         children.entry(p).or_default().push(child);
//     }
//     // Sort children for stable output.
//     for child_list in children.values_mut() {
//         child_list.sort();
//     }

//     // Build lookups for ring closures.
//     let mut open_map: BTreeMap<NodeIndex, Vec<&RingClosure>> = BTreeMap::new();
//     let mut close_map: BTreeMap<NodeIndex, Vec<&RingClosure>> = BTreeMap::new();
//     for rc in ring_closures.iter() {
//         open_map.entry(rc.opening).or_default().push(rc);
//         close_map.entry(rc.closing).or_default().push(rc);
//     }

//     // Recursively generate SMILES with a modified DFS.
//     fn dfs_tree(
//         graph: &MoleculeGraph,
//         current: NodeIndex,
//         children: &BTreeMap<NodeIndex, Vec<NodeIndex>>,
//         open_map: &BTreeMap<NodeIndex, Vec<&RingClosure>>,
//         close_map: &BTreeMap<NodeIndex, Vec<&RingClosure>>,
//         is_root: bool,
//     ) -> String {
//         // Start with the current atom's symbol.
//         let mut s = graph[current].smiles_symbol();

//         // Insert ring closure markers that open at this atom.
//         if let Some(rcs) = open_map.get(&current) {
//             for rc in rcs {
//                 let bond = match rc.bond {
//                     Bond::Single => "",
//                     Bond::Double => "=",
//                     Bond::Triple => "#",
//                     Bond::Aromatic => "",
//                 };
//                 s.push_str(bond);
//                 s.push_str(&format_ring(rc.digit));
//             }
//         }
//         // Insert ring closures that close at this atom after the atom symbol.
//         if let Some(rcs) = close_map.get(&current) {
//             for rc in rcs {
//                 let bond = match rc.bond {
//                     Bond::Single => "",
//                     Bond::Double => "=",
//                     Bond::Triple => "#",
//                     Bond::Aromatic => ":",
//                 };
//                 s.push_str(bond);
//                 s.push_str(&format_ring(rc.digit));
//             }
//         }
//         // Determine how many children this node has in the spanning tree.
//         let child_list = children.get(&current);
//         let num_children = child_list.map_or(0, |v| v.len());
//         // Compute the total degree in the spanning tree for a non-root node.
//         // For the root, we force branch notation for every child.
//         let force_branch = if is_root { true } else { num_children > 1 };

//         // Process children.
//         if let Some(child_nodes) = child_list {
//             let mut branch_str = String::new();
//             if force_branch {
//                 // For each child, output as a branch (with bond symbol) enclosed in parentheses.
//                 for &child in child_nodes.iter() {
//                     branch_str.push_str(&format!(
//                         "({}{})",
//                         bond_str(graph, current, child),
//                         dfs_tree(graph, child, children, open_map, close_map, false)
//                     ));
//                 }
//             } else {
//                 // Only one child: inline the continuation.
//                 let child = child_nodes[0];
//                 branch_str.push_str(&bond_str(graph, current, child));
//                 branch_str.push_str(&dfs_tree(
//                     graph, child, children, open_map, close_map, false,
//                 ));
//             }
//             s.push_str(&branch_str);
//         }
//         s
//     }

//     dfs_tree(graph, root, &children, &open_map, &close_map, true)
// }

// /// Returns the bond symbol for the edge between nodes `a` and `b`.
// fn bond_str(graph: &MoleculeGraph, a: NodeIndex, b: NodeIndex) -> String {
//     if let Some(edge) = graph.find_edge(a, b) {
//         let bond = graph.edge_weight(edge).unwrap();
//         match bond {
//             Bond::Single => "".to_string(),
//             Bond::Double => "=".to_string(),
//             Bond::Triple => "#".to_string(),
//             Bond::Aromatic => ":".to_string(),
//         }
//     } else {
//         "".to_string()
//     }
// }

// /// Formats a ring closure digit according to SMILES rules.
// fn format_ring(digit: usize) -> String {
//     if digit < 10 {
//         digit.to_string()
//     } else {
//         format!("%{}", digit)
//     }
// }

/// Second pass: Generate SMILES from the spanning tree, inserting ring–closure markers
/// and forcing branch notation when the current atom is a branching center.
fn generate_smiles_from_tree(
    graph: &MoleculeGraph,
    root: NodeIndex,
    parent: &BTreeMap<NodeIndex, NodeIndex>,
    ring_closures: &mut Vec<RingClosure>,
) -> String {
    // First, assign digits to ring closures.
    ring_closures.sort_by_key(|rc| (rc.opening.index(), rc.closing.index()));
    for (i, rc) in ring_closures.iter_mut().enumerate() {
        rc.digit = i + 1;
    }

    // Build the children mapping from the parent mapping.
    let mut children: BTreeMap<NodeIndex, Vec<NodeIndex>> = BTreeMap::new();
    for (&child, &p) in parent.iter() {
        children.entry(p).or_default().push(child);
    }
    // Sort children for stable output.
    for child_list in children.values_mut() {
        child_list.sort();
    }

    // Build lookups for ring closures.
    let mut open_map: BTreeMap<NodeIndex, Vec<&RingClosure>> = BTreeMap::new();
    let mut close_map: BTreeMap<NodeIndex, Vec<&RingClosure>> = BTreeMap::new();
    for rc in ring_closures.iter() {
        open_map.entry(rc.opening).or_default().push(rc);
        close_map.entry(rc.closing).or_default().push(rc);
    }

    // Recursive DFS that sorts branches using tie breakers.
    fn dfs_tree(
        graph: &MoleculeGraph,
        current: NodeIndex,
        children: &BTreeMap<NodeIndex, Vec<NodeIndex>>,
        open_map: &BTreeMap<NodeIndex, Vec<&RingClosure>>,
        close_map: &BTreeMap<NodeIndex, Vec<&RingClosure>>,
        is_root: bool,
    ) -> String {
        // Start with the current atom’s symbol.
        let mut s = graph[current].smiles_symbol();

        // Insert ring closure markers that open at this atom.
        if let Some(rcs) = open_map.get(&current) {
            for rc in rcs {
                let bond = match rc.bond {
                    Bond::Single => "",
                    Bond::Double => "=",
                    Bond::Triple => "#",
                    Bond::Aromatic => "",
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

        // Check if there are children (branches) from this node.
        if let Some(child_nodes) = children.get(&current) {
            // For each child, compute its “branch string” (using the same DFS) and also record
            // the bond that connects the parent to that child.
            let mut branches: Vec<(String, String, usize)> = Vec::new();
            for &child in child_nodes {
                let branch_smiles = dfs_tree(graph, child, children, open_map, close_map, false);
                let bond = bond_str(graph, current, child);
                let len = branch_smiles.len();
                branches.push((branch_smiles, bond, len));
            }
            // Sort the branches using tie breakers:
            //   1. Compare the branch SMILES strings,
            //   2. then the bond strings,
            //   3. then (if needed) the length of the branch string.
            branches.sort_by(|a, b| {
                let cmp = a.0.cmp(&b.0);
                if cmp != std::cmp::Ordering::Equal {
                    cmp
                } else {
                    let cmp_bond = a.1.cmp(&b.1);
                    if cmp_bond != std::cmp::Ordering::Equal {
                        cmp_bond
                    } else {
                        a.2.cmp(&b.2)
                    }
                }
            });

            // Determine if branch notation is forced.
            // For the root, we always use branches; for non-root nodes, if there is more than one child,
            // we enclose each branch in parentheses.
            let force_branch = if is_root { true } else { child_nodes.len() > 1 };

            if force_branch {
                for (branch_smiles, bond, _len) in branches {
                    s.push('(');
                    s.push_str(&bond);
                    s.push_str(&branch_smiles);
                    s.push(')');
                }
            } else {
                if let Some((branch_smiles, bond, _len)) = branches.into_iter().next() {
                    s.push_str(&bond);
                    s.push_str(&branch_smiles);
                }
            }
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
pub fn molecule_to_smiles(graph: &MoleculeGraph) -> Result<String> {
    let mut graph = graph.clone();
    graph = graph.morgan_canonize();
    // graph = graph.canonize();
    let start = graph.node_indices().next();
    match start {
        Some(start) => {
            let (parent, mut ring_closures) = compute_spanning_tree_and_ring_closures(&mut graph, start);
            Ok(generate_smiles_from_tree(&graph, start, &parent, &mut ring_closures))
        }
        None => {
            return Err(anyhow::anyhow!("Empty graph"));
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ElementType::*;

    // Assuming visualize_graph is implemented as per the previous assistant's message
    // Ensure you have the visualize_graph function in your parser.rs
    #[test]
    fn canonize() {
        let org1 = OrganicMolecule::from_smiles(&canonize_smiles("c1ccccc1R").unwrap()).unwrap();
        let org2 = OrganicMolecule::from_smiles(&canonize_smiles("c1ccc(R)cc1").unwrap()).unwrap();
        let org3 = OrganicMolecule::from_smiles(&canonize_smiles("c1cc(R)ccc1").unwrap()).unwrap();

        let smiles1 = org1.to_smiles().unwrap();
        let smiles2 = org2.to_smiles().unwrap();
        let smiles3 = org3.to_smiles().unwrap();

        println!("SMILES 1: {}", smiles1);
        println!("SMILES 2: {}", smiles2);
        println!("SMILES 3: {}", smiles3);
    }

    #[test]
    fn test_parse_ethanol() {
        let smiles = "CCO"; // Ethanol
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        // visualize_graph(&molecule, "ethanol.dot", Some("ethanol.png"))
        //     .expect("Failed to visualize graph");

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
        // visualize_graph(&molecule, "cyclohexane.dot", Some("cyclohexane.png"))
        //     .expect("Failed to visualize graph");
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
        // println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_benzene() {
        let smiles = "c1ccccc1"; // Cyclohexane
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        // visualize_graph(&molecule, "cyclohexane.dot", Some("cyclohexane.png"))
        //     .expect("Failed to visualize graph");
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
        // println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_ciprofloxacin() {
        let smiles = "C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O".to_string(); // Ciprofloxacin

        let molecule = parse_smiles(&smiles).expect("Failed to parse SMILES");
        // visualize_graph(&molecule, "ciprofloxacin.dot", Some("ciprofloxacin.png"))
        //     .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
        // println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_1_choloro_ethane() {
        let smiles = "CCCl".to_string(); // Ciprofloxacin

        let molecule = parse_smiles(&smiles).expect("Failed to parse SMILES");
        // visualize_graph(&molecule, "1-choloro-ethane.dot", Some("1-choloro-ethane.png"))
        //     .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
        // println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_isobutane() {
        let smiles = "CC(C)C"; // Isobutane
        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
        // visualize_graph(&molecule, "isobutane.dot", Some("isobutane.png"))
        //     .expect("Failed to visualize graph");

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
        // println!("Name: {}", iupac_name(&molecule));
    }

    #[test]
    fn test_parse_phenylalanine() {
        let smiles = "NC(Cc1ccccc1)C(=O)O".to_string(); // Phenylalanine

        let molecule = parse_smiles(&smiles).expect("Failed to parse SMILES");
        visualize_graph(&molecule, "phenylalanine.dot", Some("phenylalanine.png"))
            .expect("Failed to visualize graph");
        println!("Graph: {:#?}", molecule);
    }
}
