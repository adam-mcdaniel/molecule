use actix_web::{web, App, HttpResponse, HttpServer, Responder};
use anyhow::Result;
use base64::{engine::general_purpose, Engine as _};
use molecule::*; // Assumes your molecule module exports OrganicMolecule and init_logging.
use serde::Deserialize;
use tracing::*;

#[derive(Deserialize)]
struct MoleculeInput {
    input: String,
    mode: String, // "smiles", "iupac", or "sexpr"
}

const GLOBAL_STYLES: &str = r#"
/* Global styles */
body {
		font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
		background: #f5f7fa;
		color: #333;
		margin: 20px;
		padding: 0;
}
.container {
		max-width: 800px;
		margin: 50px auto;
		background: #ffffff;
		padding: 40px;
		border-radius: 8px;
		box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
}
h1 {
		font-size: 2.5rem;
		margin-bottom: 20px;
		color: #222;
		text-align: center;
}
form {
		display: flex;
		flex-direction: column;
		gap: 20px;
}
label {
		font-weight: bold;
		margin-bottom: 5px;
}
input[type="text"] {
		padding: 10px;
		font-size: 1rem;
		border: 1px solid #ccc;
		border-radius: 4px;
}
input[type="radio"] {
		margin-right: 5px;
}
input[type="submit"] {
		padding: 12px;
		background-color: #007BFF;
		color: #fff;
		border: none;
		border-radius: 4px;
		font-size: 1.1rem;
		cursor: pointer;
		transition: background-color 0.3s ease;
}
input[type="submit"]:hover {
		background-color: #0056b3;
}
img {
		display: block;
		max-width: 100%;
		height: auto;
		margin: 30px auto;
		border: 1px solid #ddd;
		border-radius: 4px;
		box-shadow: 0 2px 5px rgba(0,0,0,0.1);
}
a {
		color: #007BFF;
		text-decoration: none;
		font-weight: bold;
}
a:hover {
		text-decoration: underline;
}
.alert {
		padding: 15px;
		background-color: #f44336;
		color: white;
		margin-bottom: 20px;
		border-radius: 4px;
}
.radio-group {
	display: flex;
	flex-direction: column;
	gap: 0.5rem;
}
input[type="text"] {
    width: 100%; /* Add this line */
    padding: 10px;
    font-size: 1rem;
    border: 1px solid #ccc;
    border-radius: 4px;
}
"#;

// Some example data; replace with your real data
fn smiles_examples() -> Vec<(&'static str, &'static str)> {
    vec![
		("CCO", "Alcohol🍷"),
		("CC(=O)O", "Vinegar🍾"),
		("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Caffeine☕"),
		("NCCc1cc(O)c(O)cc1", "Dopamine🧠"),
		("C1=CC2=C(C=C1O)C(=CN2)CCN", "Serotonin🧠"),
		("c1(OC(C)(=O))(:c:c:c:c:c:1C(O)(=O))", "Aspirin💊"),
		("CC(=O)Nc1ccc(O)cc1", "Tylenol💊"),
		("C(N(C)(CCC1=CNc2:c:c:c:c:c1:2)", "DMT💊"),
		("CNC(C)Cc1ccccc1", "Methamphetamine💊"),
		("O[C@@H](c1ccccc1)[C@@H](NC)C", "Pseudoephedrine💊"),
        ("O=C(O)C(N)C", "Alanine (Amino Acid)🧬"),
		("NC(C(=O)O)CCCNC(=N)N", "Arginine (Amino Acid)🧬"),
		("C([C@@H](C(=O)O)N)C(=O)N", "Asparagine (Amino Acid)🧬"),
        ("NC(C(=O)O)CC(=O)O", "Aspartic Acid (Amino Acid)🧬"),
		("NC(C(=O)O)CS", "Cysteine (Amino Acid)🧬"),
		("C(CC(=O)N)C(C(=O)O)N", "Glutamine (Amino Acid)🧬"),
		("NC(C(=O)O)CCC(=O)O", "Glutamate (Amino Acid)🧬"),
		("NCC(=O)O", "Glycine (Amino Acid)🧬"),
		("C1=C(NC=N1)CC(C(=O)O)N", "Histidine (Amino Acid)🧬"),
		("N[C@H](C(=O)O)[C@H](CC)C", "Isoleucine (Amino Acid)🧬"),
		("NC(C(=O)O)CC(C)C", "Leucine (Amino Acid)🧬"),
		("C(CCN)CC(C(=O)O)N", "Lysine (Amino Acid)🧬"),
		("NC(C(=O)O)CCSC", "Methionine (Amino Acid)🧬"),
		("C(C(C(O)(=O))(N))(c1:c:c:c:c:c:1)", "Phenylalanine (Amino Acid)🧬"),
		("OC(=O)C1CCCN1", "Proline (Amino Acid)🧬"),
		("NC(C(=O)O)CO", "Serine (Amino Acid)🧬"),
		("NC(C(=O)O)C(C)O", "Threonine (Amino Acid)🧬"),
		("C1(=C(c2:c(:c:c:c:c:2)(N1))(CC(C(=O)(O))(N)))", "Tryptophan (Amino Acid)🧬"),
		("N[C@@H](Cc1ccc(O)cc1)C(O)=O", "Tyrosine (Amino Acid)🧬"),
		("CC(C)[C@@H](C(=O)O)N", "Valine (Amino Acid)🧬"),
    ]
}

fn iupac_examples() -> Vec<(&'static str, &'static str)> {
    vec![
		("2-methyl-1,3,5-trinitrobenzene", "TNT💣"),
        ("benzene", "Benzene🔥"),
		("cyclohexane", "Cyclohexane 🛢️"),
		("ethanol", "Ethanol🍺"),
		("methanol", "Methanol (Wood Alcohol)🔥"),
		("acetic acid", "Acetic Acid🍾"),
		("toluene", "Toluene💧"),
		("phenol", "Phenol (Anti Septic)🧪"),
		("cyclohexylbenzene", "Cyclohexylbenzene🔥"),
		("methylbenzene", "Methylbenzene🔥"),
    ]
}

fn markush_examples() -> Vec<(&'static str, &'static str)> {
    vec![
        ("benzene {CC(R)R' ammonia \"formic acid\"}", "Phenylalanine (Amino Acid)🧬"),
		("benzene benzene", "Phenylbenzene🌱"),
		("Cc1c(R)cc(R)cc1R N(=O)O N(=O)O N(=O)O", "TNT💣"),
    ]
}

/// Generate an HTML <option> list with <optgroup> for SMILES, IUPAC, and Markush.
fn generate_dropdown_options(
    smiles: &[(&str, &str)],
    iupac: &[(&str, &str)],
    markush: &[(&str, &str)],
) -> String {
    let mut html = String::new();

    // Always start with a "no selection" option
    html.push_str(r#"<option value="">--Select a preset--</option>"#);

    // SMILES group
    html.push_str(r#"<optgroup label="SMILES">"#);
    for (code, name) in smiles {
        html.push_str(&format!(
            r#"<option value="{val}" data-mode="smiles">{desc}</option>"#,
            val = code.replace("\"", "&quot;"),
            desc = name
        ));
    }
    html.push_str("</optgroup>");

    // IUPAC group
    html.push_str(r#"<optgroup label="IUPAC">"#);
    for (code, name) in iupac {
        html.push_str(&format!(
            r#"<option value="{val}" data-mode="iupac">{desc}</option>"#,
            val = code.replace("\"", "&quot;"),
            desc = name
        ));
    }
    html.push_str("</optgroup>");

    // Markush group
    html.push_str(r#"<optgroup label="Markush">"#);
    for (code, name) in markush {
        html.push_str(&format!(
            r#"<option value="{val}" data-mode="markush">{desc}</option>"#,
            val = code.replace("\"", "&quot;"),
            desc = name
        ));
    }
    html.push_str("</optgroup>");

    html
}

/// Renders the input form.
async fn index() -> impl Responder {
    // Generate the dropdown HTML from your lists
    let dropdown_html = generate_dropdown_options(
        &smiles_examples(),
        &iupac_examples(),
        &markush_examples(),
    );

    let html = format!(
        r#"
		<!DOCTYPE html>
		<html>
			<head>
				<meta charset="utf-8">
				<title>Molecule Namer + Visualizer</title>
				<style>
					{GLOBAL_STYLES}
				</style>
                <script>
                    function updateFromDropdown(selectElem) {{
                        // Get the selected option's value and data-mode attribute.
                        var selected = selectElem.options[selectElem.selectedIndex];
                        var inputField = document.getElementById('input');
                        inputField.value = selected.value;
                        
                        var mode = selected.getAttribute('data-mode');
                        if(mode) {{
                            // Check the radio button corresponding to the mode.
                            var radio = document.querySelector('input[name="mode"][value="' + mode + '"]');
                            if(radio) {{
                                radio.checked = true;
                            }}
                        }}
                    }}

                    // Optional: if you want to filter the dropdown when radio buttons change,
                    // you can add another function that hides options not matching the selected mode.
                </script>
			</head>
			<body>
				<div class="container">
					<h1>Molecule Namer + Visualizer</h1>
					<form action="/visualize" method="post">
						<label for="input">
							Enter SMILES, IUPAC, or a Markush structure:
							<input type="text" id="input" name="input" size="100" required>
						</label>

                        <!-- Dropdown to pick from predefined structures -->
                        <label for="preset">
							Or pick a preset structure:
							<select id="preset" name="preset" onchange="updateFromDropdown(this);">
								{dropdown_html}
							</select>
						</label>

						<div class="radio-group">
							<label>
								<input type="radio" name="mode" value="smiles" checked>
								SMILES
							</label>
							<label>
								<input type="radio" name="mode" value="iupac">
								IUPAC
							</label>
							<label>
								<input type="radio" name="mode" value="markush">
								Markush Editor
							</label>
						</div>
						<input type="submit" value="Name & Visualize">
					</form>
				</div>
			</body>
		</html>
		"#
    );
    HttpResponse::Ok().content_type("text/html").body(html)
}

fn error_page(error: impl std::fmt::Display) -> HttpResponse {
    let html = format!(
        r#"
				<!DOCTYPE html>
				<html>
					<head>
						<meta charset="utf-8">
						<title>Error</title>
						<style>
							{GLOBAL_STYLES}
						</style>
					</head>
					<body>
						<div class="container">
							<h1>Error</h1>
							<p>{}</p>
							<a href="/">Go back</a>
						</div>
					</body>
				</html>
		"#,
        error.to_string().replace("\n", "<br>") // Replace newlines with <br> for HTML display
    );
    HttpResponse::Ok().content_type("text/html").body(html)
}

/// Processes the form submission, converts the input to a molecule,
/// calls to_iupac() when SMILES is entered, and returns an HTML page
/// embedding the generated image and the IUPAC name.
/// On error, an alert is shown.
async fn visualize(form: web::Form<MoleculeInput>) -> impl Responder {
    // Parse the molecule based on the selected mode.
    // let molecule_result = if form.mode == "smiles" {
    //     OrganicMolecule::from_smiles(&form.input)
    // } else {
    //     OrganicMolecule::from_iupac(&form.input)
    // };
    let molecule_result = match form.mode.as_str() {
        "smiles" => Substituent::from_smiles(&form.input).map(|s| OrganicMolecule::from(s)),
        "iupac" => OrganicMolecule::from_iupac(&form.input),
        "markush" => OrganicMolecule::from_sexpr(&form.input),
        _ => Err(anyhow::anyhow!("Invalid mode")),
    };

    info!("Molecule result: {:?}", molecule_result);

    let mut molecule = match molecule_result {
        Ok(m) => m,
        Err(e) => {
            error!("Error parsing molecule: {:?}", e);
            return error_page(format!("Error parsing molecule: {:#}", e));
        }
    };

    info!("Molecule: {:#?}", molecule);
    // If the user entered SMILES, try to generate the IUPAC name.
    let iupac_display = if form.mode == "smiles" {
		// Try to lookup the supplied SMILES in the database
		match lookup_smiles_to_iupac(&form.input) {
			Some(iupac) => iupac,
			None => {
				// If not found, try to generate the IUPAC name
				match molecule.to_iupac() {
					Ok(name) => name,
					Err(e) => format!("Error converting to IUPAC: {:#}", e),
				}
			}
		}
    } else if form.mode != "iupac" {
        match molecule.to_iupac() {
            Ok(name) => name,
            Err(e) => format!("Error converting to IUPAC: {:#}", e),
        }
    } else {
        // For IUPAC input, we display the input itself.
        form.input.clone()
    };

    let smiles_display = if form.mode != "smiles" {
        match molecule.to_smiles() {
            Ok(smiles) => smiles,
            Err(e) => format!("Error converting to SMILES: {:#}", e),
        }
    } else {
        // For SMILES input, we display the input itself.
        canonize_smiles(&form.input.clone()).unwrap_or(form.input.clone())
    };

    // Generate the PNG image of the molecule.
    // We assume OrganicMolecule has a method `visualize("filename") -> Result<()>`
    if let Err(e) = molecule.visualize("output.png") {
        return error_page(format!("Error generating image: {:#}", e));
    }

    // Read the generated image.
    let png_bytes = match std::fs::read("output.png") {
        Ok(b) => b,
        Err(e) => {
            return error_page(format!("Error reading generated image: {:#}", e));
        }
    };

    // Visualize the substituents if the mode is "markush"
    // let substituent_count = if form.mode == "markush" {
	// 	molecule.visualize_substituents().unwrap_or(0)
    // } else {
	// 	0
	// };

	// Display the number of substituents
	let substituent_smiles_and_iupac = if form.mode == "markush" {
		molecule.substituents_recursive().iter().map(|s| (s.to_smiles().unwrap_or_default(), s.to_iupac().unwrap_or_default())).collect::<Vec<_>>()
	} else {
		vec![]
	};

	let substituent_png_bytes = substituent_smiles_and_iupac.iter().map(|(smiles, _)| {
		let substituent = Substituent::from_smiles(smiles).unwrap_or_else(|_| {
			error!("Error parsing substituent: {}", smiles);
			Substituent::from_smiles("R").unwrap()
		});
		if let Err(e) = substituent.visualize("substituent.png") {
			error!("Error generating image for substituent {}: {:#}", smiles, e);
			return vec![];
		}
		std::fs::read("substituent.png").unwrap_or_else(|e| {
			error!("Error reading generated image for substituent {}: {:#}", smiles, e);
			vec![]
		})
	}).collect::<Vec<_>>();
	// Encode substituent image bytes to base64.
	let substituent_encoded = substituent_png_bytes.iter().map(|bytes| general_purpose::STANDARD.encode(bytes)).collect::<Vec<_>>();

	// Now create the subsection of the page dedicated to substituents
	let substituent_html = if form.mode == "markush" {
		let mut result = String::from("<br/><h2>Substituents</h2><hr/>");
		// Create a subsubsection for each substituent
		for (i, ((smiles, iupac), encoded)) in substituent_smiles_and_iupac.iter().zip(substituent_encoded.iter()).enumerate() {
			if i == 0 {
				continue;
			}
			result.push_str(&format!(
				r#"
					<div class="subsubsection">
						<h3>Substituent #{i}</h3>
						<p><strong>SMILES:</strong> {smiles}</p>
						<p><strong>IUPAC Name:</strong> {iupac}</p>
						<img src="data:image/png;base64,{encoded}" alt="Substituent Image"/>
					</div>
				"#,
			));
		}
		result
	} else {
		String::new()
	};

    // Encode image bytes to base64.
    let encoded = general_purpose::STANDARD.encode(&png_bytes);
    let html = format!(
        r#"
		<!DOCTYPE html>
		<html>
			<head>
				<meta charset="utf-8">
				<title>Molecule Name + Visualization</title>
				<style>
					body {{
							font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
							background: #f5f7fa;
							color: #333;
							margin: 20px;
							padding: 0;
					}}
					.container {{
							max-width: 800px;
							margin: 50px auto;
							background: #ffffff;
							padding: 40px;
							border-radius: 8px;
							box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
					}}
					h1 {{
							text-align: center;
					}}
					img {{
							display: block;
							max-width: 100%;
							height: auto;
							margin: 30px auto;
							border: 1px solid #ddd;
							border-radius: 4px;
							box-shadow: 0 2px 5px rgba(0,0,0,0.1);
					}}
					a {{
							color: #007BFF;
							text-decoration: none;
							font-weight: bold;
					}}
					a:hover {{
							text-decoration: underline;
					}}
				</style>
			</head>
			<body>
				<div class="container">
					<h1>Molecule Visualization</h1>
					<p><strong>SMILES:</strong> {smiles}</p>
					<p><strong>IUPAC Name:</strong> {iupac}</p>
					<img src="data:image/png;base64,{encoded}" alt="Molecule Image"/><br>
					{substituent_html}
					<a href="/">Go back</a>
				</div>
			</body>
		</html>
		"#,
        encoded = encoded,
        iupac = iupac_display,
        smiles = smiles_display,
		substituent_html = substituent_html,
    );

    HttpResponse::Ok().content_type("text/html").body(html)
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    // Initialize logging and any molecule setup.
    init_logging("info");

    HttpServer::new(|| {
        App::new()
            .route("/", web::get().to(index))
            .route("/visualize", web::post().to(visualize))
    })
    .bind("127.0.0.1:8080")?
    .run()
    .await
}
