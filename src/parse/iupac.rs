
// // parse.rs

// use nom::{
//     branch::alt,
//     bytes::complete::{is_not, tag, take_while, take_while1},
//     character::complete::{alpha1, digit1, multispace0},
//     combinator::{all_consuming, map, map_res, opt, peek},
//     multi::{many0, many1, separated_list0},
//     sequence::{delimited, pair, preceded, separated_pair, terminated},
//     IResult,
// };

// use crate::*;

// /// This table captures a few standard Greek/Latin prefixes for chain length.
// /// You can extend it arbitrarily, or compose e.g. "undec" (11), "dodec" (12), "icos" (20), etc.
// static CHAIN_LENGTH_PREFIXES: &[(&str, usize)] = &[
//     ("meth", 1),
//     ("eth", 2),
//     ("prop", 3),
//     ("but", 4),
//     ("pent", 5),
//     ("hex", 6),
//     ("hept", 7),
//     ("oct", 8),
//     ("non", 9),
//     ("dec", 10),
//     ("undec", 11),   // example
//     ("dodec", 12),   // example
//     ("icos", 20),    // example
//     // etc. 
// ];

// /// Some ring indicators that might appear at the start of a parent name, e.g. "cyclo".
// /// You could similarly support "bicyclo", "spiro", "norborn", "azaspiro", etc.
// static RING_INDICATORS: &[&str] = &[
//     "cyclo",
//     "bicyclo",
//     "spiro",
//     // ...
// ];

// /// A minimal set of suffixes for functional groups. 
// /// In reality, you might parse them with more advanced grammar 
// /// (e.g. "carboxylic acid", "carbaldehyde", etc.).
// static FUNCTIONAL_GROUP_SUFFIXES: &[&str] = &[
//     "oic acid",
//     "al",
//     "one",
//     "ol",
//     "amine",
//     "nitrile",
//     "amide",
//     "thiol",
//     "sulfonic acid",
//     "peroxide",
//     "oate", // e.g. ester
// ];

// static FUNCTIONAL_SUFFIX_PATTERNS: &[&str] = &[
//     // store in decreasing length order so "oic acid" is matched before "acid"
//     "oic acid",
//     "carboxylic acid",
//     "al",
//     "one",
//     "ol",
//     "amine",
//     "amide",
//     "nitrile",
//     "thiol",
//     "sulfonic acid",
//     "peroxide",
//     "oate",
//     // etc.
// ];


// // -----------------------------------------------------------------------
// // Utilities / Basic Parsers
// // -----------------------------------------------------------------------

// /// Consumes any whitespace (if present).
// fn ws(input: &str) -> IResult<&str, &str> {
//     multispace0(input)
// }

// /// Parse an unsigned integer (e.g. "123") -> usize
// fn parse_usize(input: &str) -> IResult<&str, usize> {
//     map_res(digit1, |digits: &str| digits.parse::<usize>())(input)
// }

// /// Attempt to parse a ring indicator from a known table, returning (true, leftover) if recognized.
// fn parse_ring_indicator(input: &str) -> IResult<&str, bool> {
//     // We'll try each ring indicator. If none match, we return false.
//     // Because we want to be "general," we won't do alt(tag("cyclo"), tag("bicyclo"),...).
//     // Instead, we can do a "longest match" search in RING_INDICATORS.
//     // A simple approach is to check each item in RING_INDICATORS in descending length order.
//     // For demonstration, let's do a naive approach with find:
//     for indicator in RING_INDICATORS {
//         if input.starts_with(indicator) {
//             let leftover = &input[indicator.len()..];
//             return Ok((leftover, true));
//         }
//     }
//     // If none matched:
//     Ok((input, false))
// }

// /// Parse a chain length prefix by scanning the front of `input` 
// /// for any known CHAIN_LENGTH_PREFIXES. 
// /// If found, accumulate the chain length. 
// /// You can do repeated matches to handle "pentadec" => 5 + 10 = 15, etc.
// /// This is the general approach to “composing” numeric prefixes.
// fn parse_chain_length_prefix(input: &str) -> IResult<&str, usize> {
//     let mut leftover = input;
//     let mut total_length = 0;

//     // Keep matching while we find a prefix
//     'outer: loop {
//         let mut matched_something = false;
//         for &(prefix, length) in CHAIN_LENGTH_PREFIXES {
//             if leftover.starts_with(prefix) {
//                 leftover = &leftover[prefix.len()..];
//                 total_length += length;
//                 matched_something = true;
//                 break;
//             }
//         }
//         if !matched_something {
//             break 'outer;
//         }
//     }

//     if total_length == 0 {
//         // We didn't match any prefix. 
//         // Return success with 0 => means "didn't parse any chain length prefix," leftover is the same.
//         Ok((leftover, 0))
//     } else {
//         Ok((leftover, total_length))
//     }
// }

// /// Parse a single saturation or unsaturation chunk: "an", "en", "yn"
// /// For a general approach, we might parse multiple double bonds, but let's keep it simple here.
// fn parse_saturation_chunk(input: &str) -> IResult<&str, &str> {
//     alt((
//         tag("ane"), // must come before "an"
//         tag("ene"), // before "en"
//         tag("yne"), // before "yn"
//         tag("an"),
//         tag("en"),
//         tag("yn"),
//     ))(input)
// }

// /// The main "base name" parser that tries to interpret a ring indicator, 
// /// a chain-length prefix, and an optional saturation chunk. 
// /// This does *not* parse functional group suffixes. 
// /// Returns `( leftover, is_ring, chain_length, saturation )`.
// fn parse_general_parent_base(input: &str) -> IResult<&str, (bool, usize, String)> {
//     let (input, is_ring) = parse_ring_indicator(input)?;

//     // parse the numeric prefix (like "meth" => 1, "undec" => 11, etc.)
//     let (input, chain_length) = parse_chain_length_prefix(input)?;

//     // parse saturation chunk: "an", "en", or "yn". 
//     // If we fail, we'll just assume "an" by default (full saturation).
//     let (input, satur_opt) = opt(parse_saturation_chunk)(input)?;
//     let saturation = satur_opt.unwrap_or("an").to_string();

//     // Example: if we matched ring=true, length=4, satur="an", we might interpret "cyclobutan"
//     // leftover is what's not consumed. We'll pass that leftover to the suffix parser next.

//     Ok((input, (is_ring, chain_length, saturation)))
// }


// /// We'll attempt to find which suffix from FUNCTIONAL_SUFFIX_PATTERNS best matches 
// /// the start of leftover. 
// /// If none match, we return None. If one matches, we map it to a FunctionalGroup variant.
// fn parse_functional_suffix(input: &str, default_position: usize) -> IResult<&str, Option<FunctionalGroup>> {
//     // parse an optional locant first
//     let (rest, position) = opt(parse_usize)(input)?;
//     let position = position.unwrap_or(default_position);

//     // If the leftover starts with a dash, skip it
//     let (rest, _) = opt(tag("-"))(rest)?;

//     // try each known suffix
//     for pattern in FUNCTIONAL_SUFFIX_PATTERNS {
//         if rest.starts_with(pattern) {
//             let leftover = &rest[pattern.len()..];
//             // map pattern -> actual functional group
//             let fg = match *pattern {
//                 "oic acid" => FunctionalGroup::carboxylic_acid(position),
//                 "al"       => FunctionalGroup::aldehyde(position),
//                 "one"      => FunctionalGroup::ketone(position),
//                 "ol"       => FunctionalGroup::alcohol(position),
//                 "amine"    => FunctionalGroup::amine(position),
//                 "amide"    => FunctionalGroup::amide(position),
//                 "nitrile"  => FunctionalGroup::nitrile(position),
//                 "thiol"    => FunctionalGroup::thiol(position),
//                 "sulfonic acid" => FunctionalGroup::sulfonic_acid(position),
//                 "peroxide" => FunctionalGroup::peroxide(position),
//                 "oate"     => FunctionalGroup::ester(position),
//                 "carboxylic acid" => FunctionalGroup::carboxylic_acid(position),
//                 // ...
//                 _          => return Ok((leftover, None)), // fallback
//             };
//             return Ok((leftover, Some(fg)));
//         }
//     }

//     Ok((rest, None))
// }
// // -----------------------------------------------------------------------
// // Substituent Parsing
// // -----------------------------------------------------------------------


// /// A helper to parse an integer list separated by commas. E.g. "2,3,4"
// fn parse_locant_list(input: &str) -> IResult<&str, Vec<usize>> {
//     separated_list0(
//         // The separator is a comma possibly surrounded by whitespace
//         preceded(ws, tag(",")),
//         preceded(ws, parse_usize),
//     )(input)
// }

// fn parse_parent_structure_general(input: &str) -> IResult<&str, ParentStructure> {
//     // 1) parse any substituents (just as you have, with "2-chloro", "nitro", etc.)
//     let (input, subs) = parse_substituents(input)?;  // Implementation omitted for brevity

//     // 2) parse the ring indicator + chain length + saturation chunk
//     let (input, (is_ring, chain_len, satur)) = parse_general_parent_base(input)?;

//     // 3) parse functional group suffix
//     let (input, main_group) = parse_functional_suffix(input, 1)?; 
//     // position defaults to 1 unless a locant is found.

//     // 4) Build the structure: ring or chain. 
//     // For is_ring, we might do a naive ring of `chain_len` carbons, 
//     // single bonds if satur=="an", partial double if satur=="en"? 
//     // This is a simplistic approach.

//     if is_ring {
//         // e.g. build a ring with `chain_len` atoms
//         let ring = Ring {
//             atoms: (0..chain_len).map(|_| Element::C).collect(),
//             bonds: (0..chain_len).map(|_| Bond::Single).collect(),
//             is_aromatic: false,  // you'd detect aromatic differently
//         };
//         Ok((
//             input,
//             ParentStructure::Ring {
//                 ring,
//                 main_group,
//                 substituents: subs,
//                 stereocenters: vec![],
//             },
//         ))
//     } else {
//         // build a chain with `chain_len` carbons
//         Ok((
//             input,
//             ParentStructure::Chain {
//                 chain_length: if chain_len == 0 { 1 } else { chain_len }, 
//                 main_group,
//                 substituents: subs,
//                 stereocenters: vec![],
//             },
//         ))
//     }
// }

// static SUBSTITUENT_SUFFIX: &str = "yl";

// /// Example: "methyl" => chain_len=1, "ethyl" => chain_len=2, "propyl" => 3, etc.
// /// We do NOT do `alt((tag("methyl"), tag("ethyl")...))`.
// fn parse_alkyl_substituent_name(input: &str) -> IResult<&str, usize> {
//     // parse chain-length prefix
//     let (input, chain_len) = parse_chain_length_prefix(input)?;
//     // parse the final "yl"
//     let (input, _) = tag(SUBSTITUENT_SUFFIX)(input)?;
//     Ok((input, chain_len))
// }

// /// Parses multiple substituents in sequence, e.g. "2-chloro-3-methyl".
// /// Returns a Vec<Substituent>.
// pub fn parse_substituents(input: &str) -> IResult<&str, Vec<Substituent>> {
//     let mut all_subs: Vec<Substituent> = Vec::new();
//     let mut rest = input;

//     loop {
//         // Attempt to parse a single substituent (e.g. "2-chloro" or "nitro")
//         match parse_one_substituent(rest) {
//             Ok((next_input, subs)) => {
//                 // We successfully parsed one "group" of substituents, which
//                 // might include multiple locants, e.g. "2,3-dichloro" => 2 substituents
//                 all_subs.extend(subs);
//                 rest = next_input;
//                 // Some names place an extra dash between substituents,
//                 // so we consume an optional dash here before continuing.
//                 let (n_rest, _) = opt(tag("-"))(rest)?;
//                 rest = n_rest;
//             }
//             Err(_) => {
//                 // No more substituents matched; stop looping
//                 break;
//             }
//         }
//     }

//     Ok((rest, all_subs))
// }

// /// Parse a single substituent chunk:
// /// - either "2,3-chloro" => two Substituents at locant=2 and locant=3
// /// - or "nitro" => default locant=1
// fn parse_one_substituent(input: &str) -> IResult<&str, Vec<Substituent>> {
//     // Attempt path (a): parse locant-list + dash + substituent name
//     let attempt_with_locants: IResult<&str, Vec<Substituent>> = (|| {
//         let (i, locants) = parse_locant_list(input)?; // e.g. "2,3"
//         let (i, _) = tag("-")(i)?;                   // e.g. "-"
//         // parse the substituent name
//         let (i, group) = parse_substituent_name(i)?;
//         // produce a Substituent for each locant
//         let subs: Vec<Substituent> = locants.into_iter()
//             .map(|loc| Substituent {
//                 group: group.clone(),
//                 locant: loc,
//                 stereo: None,
//             })
//             .collect();
//         Ok((i, subs))
//     })();

//     match attempt_with_locants {
//         Ok(res) => Ok(res),
//         Err(_) => {
//             // Attempt path (b): no locant => parse substituent name alone => locant=1
//             let (i, group) = parse_substituent_name(input)?;
//             let sub = Substituent {
//                 group,
//                 locant: 1,
//                 stereo: None,
//             };
//             Ok((i, vec![sub]))
//         }
//     }
// }
// fn parse_substituent_name(input: &str) -> IResult<&str, FunctionalGroup> {
//     nom::branch::alt((
//         // Halogens (morphologically "chlor"+"o", "fluor"+"o", etc.)
//         map(tag("chloro"), |_| FunctionalGroup::halogen(HalogenKind::Chloro, 0)),
//         map(tag("fluoro"), |_| FunctionalGroup::halogen(HalogenKind::Fluoro, 0)),
//         map(tag("bromo"),  |_| FunctionalGroup::halogen(HalogenKind::Bromo, 0)),
//         map(tag("iodo"),   |_| FunctionalGroup::halogen(HalogenKind::Iodo, 0)),

//         // Nitro
//         map(tag("nitro"),  |_| FunctionalGroup::nitro(0)),

//         // Basic morphological parse of alkyl: parseChainPrefix + "yl" 
//         // For full generality, see prior morphological approach
//         parse_alkyl_substituent,

//         // "phenyl", "sulfanyl", "cyano", "alkoxy", etc. 
//         // If you want to parse them with a morphological approach, define a separate function. 
//         map(tag("phenyl"),   |_| FunctionalGroup::phenyl(0)),
//         map(tag("sulfanyl"), |_| FunctionalGroup::sulfanyl(0)),
//         map(tag("cyano"),    |_| FunctionalGroup::cyano(0)),
//         map(tag("alkoxy"),   |_| FunctionalGroup::ether()),

//         // fallback
//     ))(input)
// }

// /// A simple morphological approach for alkyl substituents: 
// /// E.g. "methyl" => "meth" + "yl" => chain length=1
// ///      "ethyl" => "eth" + "yl" => chain length=2
// fn parse_alkyl_substituent(input: &str) -> IResult<&str, FunctionalGroup> {
//     use nom::{sequence::tuple, bytes::complete::tag, combinator::map};

//     // We can parse the chain length prefix, then "yl"
//     let (rest, (chain_len, _)) = tuple((
//         parse_chain_length_prefix,  // e.g. "meth" => 1
//         tag("yl"),                  // final "yl"
//     ))(input)?;

//     // We can store just the name or the length. For demonstration, let's store the name as e.g. "alkyl..."
//     // Or simply store "methyl" if chain_len=1, "ethyl" if chain_len=2, etc. 
//     // Real code might do a table lookup in reverse or just keep the chain_len in a field. 
//     let name = match chain_len {
//         1 => "methyl",
//         2 => "ethyl",
//         3 => "propyl",
//         4 => "butyl",
//         5 => "pentyl",
//         6 => "hexyl",
//         7 => "heptyl",
//         8 => "octyl",
//         9 => "nonyl",
//         10 => "decyl",
//         _ => "alkyl", // fallback
//     };
//     // We'll store it in the functional group
//     Ok((rest, FunctionalGroup::Alkyl { name: name.to_string(), position: 0 }))
// }

// pub fn parse_iupac_to_parent(input: &str) -> Result<ParentStructure, String> {
//     let input = input.trim();  // optional trimming at start
//     match parse_parent_structure_general(input) {
//         Ok((rest, result)) => {
//             let leftover = rest.trim();
//             if leftover.is_empty() {
//                 // Great – we consumed the entire name
//                 Ok(result)
//             } else {
//                 // Some leftover text was never parsed
//                 Err(format!("Unparsed trailing text: '{leftover}'"))
//             }
//         }
//         Err(e) => {
//             // If the parser failed mid‐way, we convert the Nom error to a string
//             Err(format!("Parsing Error: {:?}", e))
//         }
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_parse_simple_chain() {
//         let input = "butane";
//         let result = parse_iupac_to_parent(input).expect("Should parse successfully");

//         match result {
//             ParentStructure::Chain {
//                 chain_length,
//                 main_group,
//                 substituents,
//                 ..
//             } => {
//                 assert_eq!(chain_length, 4);
//                 assert!(main_group.is_none());
//                 assert!(substituents.is_empty());
//             }
//             _ => panic!("Expected a Chain variant."),
//         }
//     }

//     #[test]
//     fn test_parse_chained_substituents() {
//         // "2-chloro-3-methylbutan-1-ol"
//         let input = "2-chloro-3-methylbutan-1-ol";
//         let result = parse_iupac_to_parent(input).expect("Should parse successfully");
//         match result {
//             ParentStructure::Chain {
//                 chain_length,
//                 main_group,
//                 substituents,
//                 ..
//             } => {
//                 // "butan" => 4 carbons
//                 assert_eq!(chain_length, 4);

//                 // The suffix is "ol" => Alcohol at position 1
//                 let mg = main_group.expect("Should have main group (alcohol).");
//                 match mg {
//                     FunctionalGroup::Alcohol { position } => {
//                         assert_eq!(position, 1);
//                     }
//                     _ => panic!("Expected an Alcohol functional group."),
//                 }

//                 // Substituents: 2-chloro, 3-methyl
//                 assert_eq!(substituents.len(), 2);

//                 // The order we get them might differ, so let's check by locant & type.
//                 // We'll just look for "chloro" at locant 2, "methyl" at locant 3.
//                 let mut found_chloro = false;
//                 let mut found_methyl = false;

//                 for sub in &substituents {
//                     match &sub.group {
//                         FunctionalGroup::Halogen { kind, position: _ } => {
//                             if *kind == HalogenKind::Chloro {
//                                 assert_eq!(sub.locant, 2);
//                                 found_chloro = true;
//                             }
//                         }
//                         FunctionalGroup::Alkyl { name, position: _ } => {
//                             if name == "methyl" {
//                                 assert_eq!(sub.locant, 3);
//                                 found_methyl = true;
//                             }
//                         }
//                         _ => {}
//                     }
//                 }

//                 assert!(found_chloro, "Should parse 2-chloro substituent.");
//                 assert!(found_methyl, "Should parse 3-methyl substituent.");
//             }
//             _ => panic!("Expected a Chain variant."),
//         }
//     }

//     #[test]
//     fn test_parse_nitrobenzene() {
//         let input = "nitrobenzene";
//         let result = parse_iupac_to_parent(input).expect("Should parse successfully");

//         match result {
//             ParentStructure::Ring {
//                 ring,
//                 main_group,
//                 substituents,
//                 ..
//             } => {
//                 // "benzene" => 6 aromatic carbons
//                 assert_eq!(ring.atoms.len(), 6);
//                 for atom in &ring.atoms {
//                     assert!(atom.is_aromatic());
//                     assert_eq!(atom.symbol(), "C");
//                 }
//                 // There's 6 aromatic bonds
//                 assert_eq!(ring.bonds.len(), 6);
//                 for b in &ring.bonds {
//                     assert_eq!(*b, Bond::Aromatic);
//                 }

//                 // No suffix group
//                 assert!(main_group.is_none());

//                 // One substituent: "nitro" at locant 1
//                 assert_eq!(substituents.len(), 1);
//                 let sub = &substituents[0];
//                 assert_eq!(sub.locant, 1);
//                 match &sub.group {
//                     FunctionalGroup::Nitro { .. } => {}
//                     _ => panic!("Expected a nitro substituent."),
//                 }
//             }
//             _ => panic!("Expected a Ring variant."),
//         }
//     }

//     #[test]
//     fn test_parse_cyclohexanone() {
//         let input = "cyclohexanone";
//         let result = parse_iupac_to_parent(input).expect("Should parse successfully");
//         println!("IUPAC: {}", result.iupac_name());

//         match result {
//             ParentStructure::Ring {
//                 ring,
//                 main_group,
//                 substituents,
//                 ..
//             } => {
//                 // "cyclohex" => 6 carbons in a ring
//                 assert_eq!(ring.atoms.len(), 6);
//                 // by default, we set them to Element::C, single bonds
//                 assert!(!ring.is_aromatic);

//                 match main_group {
//                     Some(FunctionalGroup::Ketone { position }) => {
//                         // We default to position = 1 if no locant is specified
//                         assert_eq!(position, 1);
//                     }
//                     other => panic!("Expected a ketone group, got {:?}", other),
//                 }

//                 // No substituents
//                 assert!(substituents.is_empty());
//             }
//             _ => panic!("Expected a Ring variant."),
//         }
//     }

//     #[test]
//     fn test_parse_leftover_text_error() {
//         // If there's leftover text we can't parse, we should get an Err.
//         let input = "butane leftoverstuff";
//         let result = parse_iupac_to_parent(input);
//         println!("{:?}", result);
//         assert!(result.is_err(), "Should fail due to leftover text");
//         if let Err(msg) = result {
//             assert!(
//                 msg.contains("Unparsed trailing text"),
//                 "Error should mention leftover text"
//             );
//         }
//     }
// }