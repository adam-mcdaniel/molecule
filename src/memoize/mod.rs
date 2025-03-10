use std::fs::File;
use std::io;
use crate::canonize_smiles;
use lazy_static::lazy_static;
use anyhow::Result;
use csv::{ReaderBuilder, StringRecord};
use csv::Writer;
use std::collections::{BTreeMap, BTreeSet};
use crate::Name;
use tracing::*;
use std::io::Read;

type DatabaseEntry = Name;
type CSVDatabase = Vec<(DatabaseEntry, Vec<DatabaseEntry>)>;
fn read_csv_database(csv_data: &str, smiles_column: usize, iupac_column: usize) -> CSVDatabase {
    // Create a CSV reader from the string bytes.
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(csv_data.as_bytes());
    let mut records = Vec::new();
    
    // Parse each record and push the tuple (SMILES, IUPACName) into the vector.
    for result in rdr.records() {
        let record: StringRecord = result.expect("Error reading record");
        let smiles = record.get(smiles_column).unwrap_or("").to_string();
        let iupac = record.get(iupac_column).unwrap_or("").to_string();
        if smiles.is_empty() || iupac.is_empty() {
            warn!("Skipping record with empty SMILES or IUPAC: {:?}", record);
            warn!("SMILES: '{}', IUPAC: '{}'", smiles, iupac);
            continue; // Skip records with empty SMILES or IUPAC
        }

        // Split the IUPAC name by semicolon and collect into a vector.
        let iupac: Vec<_> = iupac.split(';')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_owned())
            .map(|s| DatabaseEntry::from(s))
            .collect();

        // Push the tuple (SMILES, IUPACName) into the vector.
        records.push((smiles.into(), iupac));
    }
    records
}

lazy_static! {

    /// A static variable holding the CSV database as a vector of (SMILES, IUPACName) pairs.
    static ref EXTENDED_CSV_DATABASE: Vec<(DatabaseEntry, Vec<DatabaseEntry>)> = {
        // Include the CSV data as a string at compile time.
        // Ensure that "output.csv" is in your project's source directory or adjust the path accordingly.
        let csv_data = include_str!("extended-molecules.csv");
        
        // Read the CSV data into a string.
        read_csv_database(csv_data, 0, 1)
        // // Create a CSV reader from the string bytes.
        // let mut rdr = ReaderBuilder::new()
        //     .has_headers(true)
        //     .from_reader(csv_data.as_bytes());
        // let mut records = Vec::new();
        
        // // Parse each record and push the tuple (SMILES, IUPACName) into the vector.
        // for result in rdr.records() {
        //     let record: StringRecord = result.expect("Error reading record");
        //     let smiles = record.get(1).unwrap_or("").to_string();
        //     let iupac = record.get(3).unwrap_or("").to_string();
        //     if smiles.is_empty() || iupac.is_empty() {
        //         warn!("Skipping record with empty SMILES or IUPAC: {:?}", record);
        //         continue; // Skip records with empty SMILES or IUPAC
        //     }

        //     // Split the IUPAC name by semicolon and collect into a vector.
        //     let iupac: Vec<Name> = iupac.split(';')
        //         .map(|s| s.trim())
        //         .filter(|s| !s.is_empty())
        //         .map(|s| Name::from(adjust_iupac_name(s)))
        //         .collect();

        //     // Push the tuple (SMILES, IUPACName) into the vector.
        //     records.push((smiles.into(), iupac));
        // }
        // records
    };

    static ref LOOKUP_SMILES: BTreeMap<DatabaseEntry, DatabaseEntry> = {
        let mut map = BTreeMap::new();
        for (smiles, iupacs) in EXTENDED_CSV_DATABASE.iter() {
            // Get the shortest IUPAC name for each SMILES.
            // let iupac = iupacs.clone().into_iter()
            //     .min_by_key(|name| name.len())
            //     .unwrap_or("".to_string());
            // map.insert(smiles.clone(), iupac);
            map.insert(smiles.clone(), iupacs[0].clone());
        }
        map
    };

    static ref LOOKUP_IUPAC: BTreeMap<DatabaseEntry, DatabaseEntry> = {
        let mut map = BTreeMap::new();
        for (smiles, iupacs) in EXTENDED_CSV_DATABASE.iter() {
            for iupac in iupacs {
                // Insert the IUPAC name and corresponding SMILES into the map.
                map.insert(iupac.clone(), smiles.clone());
            }
        }
        map
    };
}

/// Adjusts a user-input IUPAC name by trimming whitespace,
/// collapsing multiple spaces to one, and removing hyphens adjacent to spaces.
fn adjust_iupac_name(input: &str) -> String {
    // Trim leading and trailing whitespace.
    let trimmed = input.trim().to_ascii_lowercase();
    // Collapse multiple whitespace characters into a single space.
    let collapsed = trimmed.split_whitespace().collect::<Vec<&str>>().join(" ");
    // Remove hyphens that appear adjacent to spaces.
    let no_adjacent_hyphens = collapsed.replace("- ", " ").replace(" -", " ");
    // Also remove hyphens at the start or end.
    no_adjacent_hyphens.trim_matches('-').to_string()
}

/// Attempts to look up the SMILES string for the given IUPAC name.
/// If a direct match is not found in the global LOOKUP_IUPAC map,
/// the function will adjust the input (normalize whitespace and hyphens)
/// and try again.
pub fn lookup_iupac_to_smiles(input: &str) -> Option<String> {
    // First, try a direct lookup.
    let iupac = DatabaseEntry::from(input.to_owned());
    if let Some(smiles) = LOOKUP_IUPAC.get(&iupac) {
        return Some(smiles.to_string());
    } else {
        // Debug print to show the input (optional)
        warn!("IUPAC input not found directly: '{}'", input);
    }
    // Try adjusting the input.
    let adjusted = adjust_iupac_name(input);
    // Debug print to show the adjusted name (optional)
    info!("Adjusted IUPAC input: '{}'", adjusted);
    let adjusted = DatabaseEntry::from(adjusted.to_owned());
    LOOKUP_IUPAC.get(&adjusted).map(|smiles| {
        info!("Found SMILES for adjusted IUPAC: '{}'", smiles);

        smiles.to_string()
    })
}


pub fn lookup_smiles_to_iupac(input: &str) -> Option<String> {
    // First, try a direct lookup.
    let smiles = DatabaseEntry::from(input.to_owned());
    if let Some(iupac) = LOOKUP_SMILES.get(&smiles) {
        info!("Found IUPAC for SMILES: '{}'", smiles);
        return Some(iupac.to_string());
    } else {
        // Debug print to show the input (optional)
        error!("SMILES input not found directly: '{}'", input);
        None
    }
}

/// Extends the CSV database by computing a canonical SMILES for each entry.
///
/// Returns a vector of tuples:
/// (canonical_smiles, original_smiles, iupac_name)
fn extend_csv_database_with_canonical(csv_data: &str) -> Result<Vec<(String, String, String)>> {
    let mut extended = Vec::new();
    let database = read_csv_database(csv_data, 1, 2);
    for (orig_smiles, iupac) in database.iter() {
        // Compute the canonical SMILES using our custom function.
        let canon_smiles = canonize_smiles(&canonize_smiles(&canonize_smiles(orig_smiles.as_ref())?)?)?;
        let iupacs: Vec<&str> = iupac.iter().map(|name| name.as_ref()).collect();
        let iupac = iupacs.join("; ");
        extended.push((canon_smiles, orig_smiles.to_string(), iupac));
    }
    Ok(extended)
}

/// Extends the CSV database by computing a canonical SMILES for each entry,
/// and then writes the extended database to a CSV file.
/// The output CSV will have three columns:
/// "Canonical_SMILES", "Original_SMILES", and "IUPAC_Name".
pub fn extend_and_write_csv_database(csv_file: &str, output_filename: &str) -> Result<()> {
    // Create an output file.
    let file = File::create(output_filename)?;
    let mut wtr = Writer::from_writer(file);

    // Read the CSV data into a string.
    let mut csv = File::open(csv_file)?;
    let mut csv_data = String::new();
    csv.read_to_string(&mut csv_data)?;

    // Write the header row.
    wtr.write_record(&["Canonical_SMILES", "IUPAC_Name"])?;

    let extended = extend_csv_database_with_canonical(&csv_data)?;
    // Write each record to the CSV.
    for (canon_smiles, _orig_smiles, iupac) in extended {
        wtr.write_record(&[canon_smiles, iupac])?;
    }
    // Ensure all data is flushed to disk.
    wtr.flush()?;
    info!("Extended CSV database written to {}", output_filename);
    Ok(())
}
