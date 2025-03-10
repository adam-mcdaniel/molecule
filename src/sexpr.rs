///! A tool for parsing S-Expressions

use nom::{
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{char, multispace0, multispace1},
    combinator::{map, cut, all_consuming, opt},
    error::{convert_error, VerboseError},
    multi::{separated_list0, many0, many1},
    sequence::{delimited, pair, preceded, terminated, tuple},
    IResult,
};

use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::str::FromStr;

// ---------------------------------------------------------------------
// Error and Location
// ---------------------------------------------------------------------

/// The error type we will use (as requested).
pub type Error<'a> = VerboseError<&'a str>;

/// A convenient alias for our IResult with that error type.
pub type Res<'a, T> = IResult<&'a str, T, Error<'a>>;

/// We will store a `Location` with each node so that we can indicate
/// where (in the original input) that node was parsed from.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Location<'a> {
    /// This will store the original input slice of this node.
    pub fragment: &'a str,
}

/// A helper that returns a `Location` for the current input position.
fn current_location(input: &str) -> Location {
    Location {
        fragment: input,
    }
}

impl Debug for Location<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "Location({} to the end)", self.fragment.len())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SExpr {
    Atom(String),
    List(Vec<SExpr>),
}

impl SExpr {
    pub fn is_list(&self) -> bool {
        matches!(self, SExpr::List(_))
    }

    pub fn is_atom(&self) -> bool {
        matches!(self, SExpr::Atom(_))
    }

    pub fn as_list(&self) -> &Vec<SExpr> {
        match self {
            SExpr::List(l) => l,
            _ => panic!("{} is not a list", self),
        }
    }

    pub fn as_atom(&self) -> &String {
        match self {
            SExpr::Atom(a) => a,
            _ => panic!("{} is not an atom", self),
        }
    }
}

// impl Debug for SExpr {
//     fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
//         match self {
//             SExpr::Atom(s) => write!(f, "{}", s),
//             SExpr::List(l) => {
//                 write!(f, "(")?;
//                 for (i, e) in l.iter().enumerate() {
//                     if i > 0 {
//                         write!(f, " ")?;
//                     }
//                     write!(f, "{:?}", e)?;
//                 }
//                 write!(f, ")")
//             }
//         }
//     }
// }

impl Display for SExpr {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            SExpr::Atom(s) => write!(f, "{}", s),
            SExpr::List(l) => {
                write!(f, "(")?;
                for (i, e) in l.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{}", e)?;
                }
                write!(f, ")")
            }
        }
    }
}

impl FromStr for SExpr {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse(s).map_err(|e| e.to_string())
    }
}


// ---------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------

/// Parses a string literal, e.g. "hello (world)"
fn parse_string_literal(input: &str) -> Res<SExpr> {
    // Allow any characters until the closing quote.
    let (input, _) = multispace0(input)?;
    let (input, content) = delimited(
        char('"'),
        // We use take_while; you might extend this with escapes.
        take_while(|c| c != '"'),
        char('"')
    )(input)?;
    Ok((input, SExpr::Atom(content.to_string())))
}

fn parse_atom(input: &str) -> Res<SExpr> {
    let (input, _) = multispace0(input)?;
    let (input, atom) = alt((
        parse_string_literal,
        // We use take_while; you might extend this with escapes.
        map(take_while(|c: char| c != '{' && c != '}' && !c.is_whitespace()), |atom: &str| {
            SExpr::Atom(atom.to_string())
        }),
    ))(input)?;
    // let (input, atom) = take_while1(|c: char| c != '{' && c != '}' && !c.is_whitespace())(input)?;
    Ok((input, SExpr::Atom(atom.to_string())))
}


fn parse_list(input: &str) -> Res<SExpr> {
    let (input, _) = multispace0(input)?;
    let (input, list) = delimited(
        char('{'),
        separated_list0(multispace1, parse_sexpr),
        preceded(multispace0, char('}')),
    )(input)?;
    Ok((input, SExpr::List(list)))
}

fn parse_sexpr(input: &str) -> Res<SExpr> {
    let (input, _) = multispace0(input)?;
    alt((parse_list, parse_atom))(input)
}

fn parse(input: &str) -> Result<SExpr, String> {
    match all_consuming(parse_sexpr)(input.trim()) {
        Ok((_, sexpr)) => Ok(sexpr),
        Err(e) => match e {
            nom::Err::Error(e) | nom::Err::Failure(e) => Err(convert_error(input, e)),
            nom::Err::Incomplete(_) => Err("incomplete".to_string()),
        }
    }
}


#[cfg(test)]
mod tests {
    #[test]
    fn test_parse() {
        assert_eq!(
            super::parse("{a b c}"),
            Ok(super::SExpr::List(vec![
                super::SExpr::Atom("a".to_string()),
                super::SExpr::Atom("b".to_string()),
                super::SExpr::Atom("c".to_string())
            ]))
        );
    }

    #[test]
    fn test_parse_nested() {
        assert_eq!(
            super::parse("{a {b c} d}"),
            Ok(super::SExpr::List(vec![
                super::SExpr::Atom("a".to_string()),
                super::SExpr::List(vec![
                    super::SExpr::Atom("b".to_string()),
                    super::SExpr::Atom("c".to_string())
                ]),
                super::SExpr::Atom("d".to_string())
            ]))
        );
    }

    #[test]
    fn test_parse_smiles() {
        assert_eq!(
            super::parse("
            C1=RC=CC=C1
            "),
            Ok(super::SExpr::Atom("C1=RC=CC=C1".to_string()))
        );
    }
    
}