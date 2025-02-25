use molecule::{parse_smiles, visualize_graph, iupac_name};

fn main() {
    let smiles = "CCO"; // Ethanol
    
    let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
    
    assert_eq!(molecule.node_count(), 3); // 2 Carbons, 1 Oxygen
    
    visualize_graph(&molecule, "ethanol.dot", Some("ethanol.png")).expect("Failed to visualize graph");

    println!("IUPAC name: {}", iupac_name(&molecule));
}