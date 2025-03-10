use molecule::*;
use anyhow::Result;

fn main() -> Result<()> {
    // init_logging("trace");


    let sexpr_molecule = OrganicMolecule::from_sexpr(r#"
        "c1ccc(R)cc1" ("R''CC(R)R'" RN {"RC(O)=O")
    "#).expect("Failed to parse S-expression");
    sexpr_molecule.visualize("sexpr_molecule.png")?;

    // "C1(=C(R)C=C(C=C1)" ("R''CC(R)R'" RN "C(O)=O")

    // let input_csv = "molecules.csv";
    // let output_csv = "extended-molecules.csv";
    // extend_and_write_csv_database(input_csv, output_csv)?;

    let phenylalanine = OrganicMolecule::from_iupac("2-amino-3-phenylpropanoic acid").expect("Failed to parse IUPAC name");
    println!("Phenylalanine: {:?}", phenylalanine);
    println!("Phenylalanine SMILES: {}", phenylalanine.to_smiles().expect("Failed to convert to SMILES"));

    let cyclohexylbenzene = OrganicMolecule::from_iupac("cyclohexylbenzene").expect("Failed to parse IUPAC name");
    println!("Cyclohexylbenzene: {:?}", cyclohexylbenzene);

    let methylbenzene = OrganicMolecule::from_iupac("methylbenzene").expect("Failed to parse IUPAC name");
    println!("Methylbenzene: {:?}", methylbenzene);

    // Visualize the molecule
    phenylalanine.visualize("phenylalanine.png")?;
    
    cyclohexylbenzene.visualize("cyclohexylbenzene.png")?;

    methylbenzene.visualize("methylbenzene.png")?;

    let strychnine = OrganicMolecule::from_smiles("O=C7N2c1ccccc1[C@@]64[C@@H]2[C@@H]3[C@@H](OC/C=C5\\[C@@H]3C[C@@H]6N(CC4)C5)C7").expect("Failed to parse SMILES");
    println!("Strychnine: {:?}", strychnine);
    // Visualize the molecule
    strychnine.visualize("strychnine.png")?;
    // Determine the IUPAC name of the molecule
    let iupac_name = phenylalanine.to_iupac().expect("Failed to convert to IUPAC");
    Ok(())
}


    /*
    let adams_custom_substituent =
        Substituent::parse_smiles("C1=C(CCCR)C=C(CCR')C=C(CR'')1").unwrap();
    println!("Adam's custom substituent: {:?}", adams_custom_substituent);

    let molecule = OrganicMolecule::from(adams_custom_substituent)
        .bind([
            // Bind a nitrogen to R
            Element::N.into(),
            // Bind a phenyl group to R'
            Substituent::phenyl_group().into(),
            // Bind another, nested bound substituent to R''
            OrganicMolecule::from(Substituent::parse_smiles("RNNNR'").unwrap())
                // Take 3 nitrogens, bind a chlorine to the first one,
                // then bind the other side to my custom substituent
                .bind([Substituent::halide(Element::Cl).into()])
                .unwrap(),
        ])
        .unwrap();

    molecule.visualize("binding-to-sites.png").unwrap();
     */

    /*
    let ethyl_pryidine_3_carboxylate =
        OrganicMolecule::parse_smiles("C(C)OC(=O)C=1C=NC=CC1").unwrap();

    // Visualize the molecule
    ethyl_pryidine_3_carboxylate
        .visualize("ethyl-pyridine-3-carboxylate.png")
        .unwrap();

    // Determine the IUPAC name of the molecule
    let iupac_name = ethyl_pryidine_3_carboxylate.to_iupac().unwrap();
    println!("IUPAC name: {}", iupac_name);

    // Determine the molecular formula of the molecule
    // println!("{}", OrganicMolecule::parse_smiles("CCOC(=O)CC(=O)OCC").unwrap().to_iupac().unwrap());
    println!(
        "{}",
        OrganicMolecule::parse_smiles("O=C(OC)C(=O)OC")
            .unwrap()
            .to_iupac()
            .unwrap()
    );
     */