use molecule::{parse_smiles, visualize_graph, OrganicMolecule};

fn main() {
    // let smiles = "COC(C)=Oasdf";
    
    // let molecule = parse_smiles(smiles).expect("Failed to parse SMILES");
    
    // // assert_eq!(molecule.node_count(), 3); // 2 Carbons, 1 Oxygen
    
    // println!("{:#?}", molecule);

    // visualize_graph(&molecule, "benzene.dot", Some("benzene.png")).expect("Failed to visualize graph");

    let benzyl_acetate = OrganicMolecule::parse_smiles("c1ccccc1COC(=O)C").unwrap();
    println!("{:#?}", benzyl_acetate);
    benzyl_acetate.visualize("benzyl_acetate.png").unwrap();

    let isoamyl_acetate = OrganicMolecule::parse_smiles("O=C(OCCC(C)C)C").unwrap();
    println!("{:#?}", isoamyl_acetate);
    isoamyl_acetate.visualize("isoamyl_acetate.png").unwrap();

    let isobutyl_formate = OrganicMolecule::parse_smiles("O=COCC(C)C").unwrap();
    println!("{:#?}", isobutyl_formate);
    isobutyl_formate.visualize("isobutyl_formate.png").unwrap();
    // println!("IUPAC name: {}", iupac_name(&molecule));
}