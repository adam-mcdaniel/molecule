// src/parser.rs

use crate::{Bond, Bond::*, Element, Element::*};
use crate::MoleculeGraph;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use std::io::Write;
use std::fmt::Write as FmtWrite;
use petgraph::prelude::EdgeRef;
use super::*;
use std::collections::HashSet;
/// Visualizes the MoleculeGraph by exporting it to a DOT format and optionally rendering it as an image.
///
/// # Arguments
///
/// * `graph` - The MoleculeGraph to visualize.
/// * `output_dot` - The path to save the DOT file.
/// * `output_image` - Optional path to save the rendered image (e.g., "molecule.png").
///
/// # Returns
///
/// * `Result<(), String>` - Ok if successful, Err with message otherwise.
pub fn visualize_graph(
    graph: &MoleculeGraph,
    output_dot: &str,
    output_image: Option<&str>,
) -> Result<(), String> {
    // Generate the DOT string with custom edge labels and styles
    let dot_string = generate_dot(graph);

    // Write the DOT string to a file
    let mut file = std::fs::File::create(output_dot)
        .map_err(|e| format!("Failed to create DOT file: {}", e))?;
    file.write_all(dot_string.as_bytes())
        .map_err(|e| format!("Failed to write to DOT file: {}", e))?;

    println!("DOT file saved to {}", output_dot);

    // If an image output path is provided, render the DOT to an image using Graphviz
    if let Some(image_path) = output_image {
        // Execute the 'dot' command to generate the image
        // Ensure Graphviz is installed and 'dot' is in the system PATH
        let status = std::process::Command::new("dot")
            .args(&["-Tpng", output_dot, "-o", image_path])
            .status()
            .map_err(|e| format!("Failed to execute Graphviz 'dot' command: {}", e))?;

        if !status.success() {
            return Err(format!(
                "Graphviz 'dot' command failed with status: {}",
                status
            ));
        }

        println!("Image rendered to {}", image_path);
    }

    Ok(())
}

/// Generates a DOT representation of the MoleculeGraph with customized node and edge attributes.
///
/// # Arguments
///
/// * `graph` - The MoleculeGraph to convert.
///
/// # Returns
///
/// * `String` - The DOT format string.
fn generate_dot(graph: &MoleculeGraph) -> String {
    // Use a String to accumulate the DOT representation.
    let mut dot_output = String::new();
    writeln!(dot_output, "graph Molecule {{").unwrap();
    writeln!(dot_output, "    layout=neato; rankdir=LR;").unwrap();
    // Allow multiple edges between the same nodes.
    writeln!(dot_output, "    multiedge=true;").unwrap();

    // Define node labels (element symbols) and colors.
    for node in graph.node_indices() {
        let element = &graph[node];
        let node_color = element_to_color(element);
        writeln!(
            dot_output,
            "    {} [label=\"{}\", fontcolor=white, shape=circle, style=filled, fillcolor={color}];",
            node.index(),
            element_to_string(element),
            color = node_color
        )
        .unwrap();
    }

    // Define edges with bond types.
    // For undirected graphs, we output each edge once (when source < target)
    // and output multiple edge statements for double or triple bonds.
    for edge in graph.edge_references() {
        let source = edge.source().index();
        let target = edge.target().index();
        let bond = edge.weight();
        let (style, penwidth, extra) = bond_to_style(bond);

        // Determine how many parallel edges to draw:
        let count = match bond {
            Bond::Double => 2,
            Bond::Triple => 3,
            _ => 1,
        };

        for _ in 0..count {
            writeln!(
                dot_output,
                "    {} -- {} [style={}, penwidth={}, {}];",
                source, target, style, penwidth, extra
            )
            .unwrap();
        }
    }

    writeln!(dot_output, "}}").unwrap();

    dot_output
}

/// Converts an Element enum to its string representation.
fn element_to_string(element: &Element) -> String {
    use Element::*;
    match element {
        C => "C".to_string(),
        CAromatic => "C".to_string(),
        H => "H".to_string(),
        O => "O".to_string(),
        OAromatic => "O".to_string(),
        N => "N".to_string(),
        NAromatic => "N".to_string(),
        Cl => "Cl".to_string(),
        Br => "Br".to_string(),
        F => "F".to_string(),
        S => "S".to_string(),
        SAromatic => "S".to_string(),
        // Add other elements as needed
    }
}

/// Assigns colors to elements for visualization.
fn element_to_color(element: &Element) -> &'static str {
    match element {
        C => "black",
        CAromatic => "black",
        H => "gray",
        O => "red",
        OAromatic => "red",
        N => "blue",
        NAromatic => "blue",
        Cl => "darkgreen",
        Br => "brown",
        F => "pink",
        S => "yellow",
        SAromatic => "yellow",
        // Add other elements with distinct colors
    }
}

/// Maps Bond types to Graphviz edge styles, pen widths, and additional attributes.
///
/// # Arguments
///
/// * `bond` - The bond type to map.
///
/// # Returns
///
/// * `(String, f64, String)` - A tuple containing the Graphviz style, penwidth, and any extra attributes.
fn bond_to_style(bond: &Bond) -> (String, f64, String) {
    match bond {
        Bond::Single => ("solid".to_string(), 2.0, "".to_string()),
        Bond::Double => ("solid".to_string(), 2.0, "".to_string()),
        Bond::Triple => ("solid".to_string(), 2.0, "".to_string()),
        Bond::Aromatic => ("dashed".to_string(), 2.0, "color=purple".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::prelude::*;

    #[test]
    fn test_element_to_string() {
        assert_eq!(element_to_string(&Element::C), "C");
        assert_eq!(element_to_string(&Element::H), "H");
        assert_eq!(element_to_string(&Element::O), "O");
        assert_eq!(element_to_string(&Element::N), "N");
        assert_eq!(element_to_string(&Element::Cl), "Cl");
        assert_eq!(element_to_string(&Element::Br), "Br");
        assert_eq!(element_to_string(&Element::F), "F");
        assert_eq!(element_to_string(&Element::CAromatic), "C");
        assert_eq!(element_to_string(&Element::OAromatic), "O");
        assert_eq!(element_to_string(&Element::NAromatic), "N");
        assert_eq!(element_to_string(&Element::SAromatic), "S");
    }

    #[test]
    fn draw_ethanol() {
        let smiles = "CCO"; // Ethanol

        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");

        assert_eq!(molecule.node_count(), 3); // 2 Carbons, 1 Oxygen

        visualize_graph(&molecule, "ethanol.dot", Some("ethanol.png")).expect("Failed to visualize graph");
    }

    #[test]
    fn draw_methyl_ethanoate() {
        let smiles = "COC(C)=O"; // Methyl ethanoate

        let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");

        assert_eq!(molecule.node_count(), 5); // 3 Carbons, 1 Oxygen, 1 Hydrogen

        visualize_graph(&molecule, "methyl_ethanoate.dot", Some("methyl_ethanoate.png")).expect("Failed to visualize graph");
    }
}