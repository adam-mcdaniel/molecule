import csv
import requests
from rdkit import Chem
import logging

# Configure logging: both console and file.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create handlers.
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)  # Adjust level as needed.

file_handler = logging.FileHandler("script.log")
file_handler.setLevel(logging.DEBUG)  # Log more detail to file.

# Create formatters and add them to handlers.
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

# Add handlers to logger.
logger.addHandler(console_handler)
logger.addHandler(file_handler)


ALLOWED_ELEMENTS = {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"}

def is_allowed_molecule(smiles):
    if '%' in smiles:
        logger.info(f"Skipping {smiles} due to massive rings.")
        return False
    
    """Return True if the molecule (from SMILES) contains only allowed atoms,
    has no isotopes, and is a single connected fragment."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Skipping {smiles} because it could not be parsed.")
        return False
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            logger.info(f"Skipping {smiles} due to isotope: {atom.GetIsotope()}")
            return False
        if atom.GetSymbol() not in ALLOWED_ELEMENTS:
            logger.info(f"Skipping {smiles} due to disallowed element: {atom.GetSymbol()}")
            return False
    frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    if len(frags) > 1:
        logger.info(f"Skipping {smiles} due to multiple fragments.")
        return False
    return True

def fetch_iupac_name_from_pubchem(cid):
    """
    Given a SMILES string, use the PubChem PUG REST API to get the IUPAC name.
    If the API call fails or no IUPAC name is found, return "Not found".
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/CSV"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            lines = response.text.splitlines()
            # CSV header assumed like "CID,IUPACName"
            if len(lines) > 1:
                parts = lines[1].split(',')
                if len(parts) >= 2:
                    return parts[1].strip()
    except Exception as e:
        logger.error(f"Error fetching IUPAC name for CID {cid}: {e}")
    return "Not found"

def main():
    sdf_file = "Compound_000000001_000500000.sdf"  # Replace with the path to your SDF file
    output_csv = "molecules.csv"
    
    # Load molecules from the SDF file.
    suppl = Chem.SDMolSupplier(sdf_file)
    compounds = []  # Will hold tuples: (CID, SMILES, IUPAC)
    count = 0
    
    for mol in suppl:
        if mol is None:
            continue

        # Retrieve the CID; skip molecule if not available.
        cid = None
        if mol.HasProp("PUBCHEM_COMPOUND_CID"):
            try:
                cid = int(mol.GetProp("PUBCHEM_COMPOUND_CID"))
            except Exception as e:
                pass
        else:
            pass
        

        if mol.HasProp("PUBCHEM_SMILES"):
            smiles = mol.GetProp("PUBCHEM_SMILES")
        else:
            smiles = Chem.MolToSmiles(mol)

        if not is_allowed_molecule(smiles):
            continue

        iupac_synonyms = []

        if mol.HasProp("PUBCHEM_IUPAC_NAME"):
            iupac_synonyms.append(mol.GetProp("PUBCHEM_IUPAC_NAME"))
        
        if mol.HasProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME"):
            iupac_synonyms.append(mol.GetProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME"))
        else:
            if cid is None:
                logger.info(f"Skipping {smiles} due to missing CID.")
                continue
            logger.info(f"Fetching IUPAC name from PubChem for {smiles}")
            iupac = fetch_iupac_name_from_pubchem(cid)
            if iupac == "Not found":
                logger.info(f"Skipping {smiles} due to missing IUPAC name.")
                continue
            else:
                iupac_synonyms.append(iupac)
                
        if mol.HasProp("PUBCHEM_IUPAC_TRADITIONAL_NAME"):
            iupac_synonyms.append(mol.GetProp("PUBCHEM_IUPAC_TRADITIONAL_NAME"))
        if mol.HasProp("PUBCHEM_IUPAC_CAS_NAME"):
            iupac_synonyms.append(mol.GetProp("PUBCHEM_IUPAC_CAS_NAME"))
        if mol.HasProp("PUBCHEM_IUPAC_OPENEYE_NAME"):
            iupac_synonyms.append(mol.GetProp("PUBCHEM_IUPAC_OPENEYE_NAME"))
        
        # Remove repeated names and join them with a semicolon.
        iupac = "; ".join(sorted(set(iupac_synonyms)))
        logger.info(f"Processing CID: {cid}, SMILES: {smiles}, IUPAC: {iupac}")

        compounds.append((cid, smiles, iupac))
        count += 1
        
        if count % 100 == 0:
            logger.info(f"Processed {count} compounds.")
            
        if count >= 50000:
            break

    # Sort compounds by CID (ascending)
    compounds.sort(key=lambda tup: tup[0])
    
    # Write results to CSV.
    with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["CID", "SMILES", "IUPAC"])
        for cid, smiles, iupac in compounds:
            writer.writerow([cid, smiles, iupac])
    
    logger.info(f"CSV file '{output_csv}' has been created with {len(compounds)} compounds.")

# ALLOWED_ELEMENTS = {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"}

# def is_allowed_molecule(smiles):
#     if '%' in smiles:
#         print(f"Skipping {smiles} due to massive rings.")
#         return False
    
#     """Return True if the molecule (from SMILES) contains only CHNOPS atoms,
#     has no isotopes, and is a single connected fragment."""
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         return False
#     # Exclude atoms with isotope information and check allowed elements.
#     for atom in mol.GetAtoms():
#         if atom.GetIsotope() != 0:
#             print(f"Skipping {smiles} due to isotope: {atom.GetIsotope()}")
#             return False
#         if atom.GetSymbol() not in ALLOWED_ELEMENTS:
#             print(f"Skipping {smiles} due to disallowed element: {atom.GetSymbol()}")
#             return False
#     # Exclude disconnected molecules.
#     frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
#     if len(frags) > 1:
#         print(f"Skipping {smiles} due to multiple fragments.")
#         return False
#     return True

# def fetch_iupac_name_from_pubchem(smiles):
#     """
#     Given a SMILES string, use the PubChem PUG REST API to get the IUPAC name.
#     If the API call fails or no IUPAC name is found, return "Not found".
#     """
#     url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/CSV"
#     try:
#         response = requests.get(url, timeout=10)
#         if response.status_code == 200:
#             lines = response.text.splitlines()
#             # CSV header assumed like "CID,IUPACName"
#             if len(lines) > 1:
#                 parts = lines[1].split(',')
#                 if len(parts) >= 2:
#                     return parts[1].strip()
#     except Exception as e:
#         print(f"Error fetching IUPAC name for {smiles}: {e}")
#     return "Not found"

# def main():
#     sdf_file = "Compound_000000001_000500000.sdf"  # Replace with the path to your SDF file
#     output_csv = "output.csv"
    
#     # Load molecules from the SDF file
#     suppl = Chem.SDMolSupplier(sdf_file)
#     compounds = []
#     count = 0
    
#     for mol in suppl:
#         if mol is None:
#             continue

#         # Retrieve the CID; skip the molecule if not available.
#         if mol.HasProp("PUBCHEM_COMPOUND_CID"):
#             try:
#                 cid = int(mol.GetProp("PUBCHEM_COMPOUND_CID"))
#             except Exception as e:
#                 pass
#         else:
#             pass
        
#         # Get the SMILES string: first try to retrieve the PubChem property, otherwise compute it.
#         if mol.HasProp("PUBCHEM_SMILES"):
#             smiles = mol.GetProp("PUBCHEM_SMILES")
#         else:
#             smiles = Chem.MolToSmiles(mol)

#         # Only include molecules with exclusively CHNOPS elements.
#         if not is_allowed_molecule(smiles):
#             continue

#         # Get the IUPAC name: try to retrieve the property, otherwise fetch from PubChem API.
#         if mol.HasProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME"):
#             # print(f"Using IUPAC name from SDF for {smiles}")
#             iupac = mol.GetProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME")
#         else:
#             # print(f"Fetching IUPAC name from PubChem for {smiles}")
#             iupac = fetch_iupac_name_from_pubchem(smiles)
#             if iupac == "Not found":
#                 print(f"Skipping {smiles} due to missing IUPAC name.")
#                 continue
#         print(f"Processing CID: {cid}, SMILES: {smiles}, IUPAC: {iupac}")
        
#         compounds.append((smiles, iupac))
#         count += 1
        
#         if count % 100 == 0:
#             print(f"Processed {count} compounds.")
            
#         if count >= 50000:
#             break

#     # Write results to CSV
#     with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(["SMILES", "IUPAC"])
#         for smiles, iupac in compounds:
#             writer.writerow([smiles, iupac])
    
#     print(f"CSV file '{output_csv}' has been created with {len(compounds)} compounds.")

if __name__ == "__main__":
    main()