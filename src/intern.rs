use std::{
    collections::HashMap,
    fmt::{Display, Debug, Formatter, Result as FmtResult},
    sync::{Arc, RwLock},
};

use lazy_static::lazy_static;

lazy_static! {
    static ref INTERNED_NAMES: RwLock<HashMap<String, Name>> = RwLock::new(HashMap::new());
}

/// A symbol that uses string interning
#[allow(clippy::derived_hash_with_manual_eq, clippy::derive_ord_xor_partial_ord)]
#[derive(Clone, Hash, Eq, Ord)]
pub struct Name(Arc<String>);

impl Name {
    /// Create a new symbol
    pub fn new(name: &str) -> Self {
        let mut symbols = INTERNED_NAMES.write().unwrap();
        if let Some(symbol) = symbols.get(name) {
            return symbol.clone();
        }

        let symbol = Name(Arc::new(name.to_string()));
        symbols.insert(name.to_string(), symbol.clone());
        symbol
    }

    /// Get the name of the symbol
    pub fn name(&self) -> &str {
        &self.0
    }

    /// Get an iterator over all symbols
    pub fn all_symbols() -> Vec<Name> {
        INTERNED_NAMES.read().unwrap().values().cloned().collect()
    }
}

impl From<&str> for Name {
    fn from(s: &str) -> Self {
        Name::new(s)
    }
}

impl From<String> for Name {
    fn from(s: String) -> Self {
        Name::new(&s)
    }
}

impl PartialEq for Name {
    fn eq(&self, other: &Self) -> bool {
        if Arc::ptr_eq(&self.0, &other.0) {
            return true;
        }
        self.0 == other.0
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for Name {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if Arc::ptr_eq(&self.0, &other.0) {
            return Some(std::cmp::Ordering::Equal);
        }
        self.0.partial_cmp(&other.0)
    }
}

impl Debug for Name {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{}", self.0)
    }
}

impl Display for Name {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "{}", self.0)
    }
}

impl AsRef<str> for Name {
    fn as_ref(&self) -> &str {
        &self.0
    }
}