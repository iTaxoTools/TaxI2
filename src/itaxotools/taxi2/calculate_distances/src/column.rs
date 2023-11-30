//! Wrapper for pandas column

use pyo3::{
    types::{IntoPyDict, PyString},
    PyAny, PyResult, Python,
};

/// Represents a column of strings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Column<'s> {
    pub strings: Box<[&'s str]>,
}

impl<'s> Column<'s> {
    /// Make a [Column] from a pandas string Series.
    pub fn new(py: Python, column: &'s PyAny) -> PyResult<Self> {
        let str_accessor = column.getattr("str")?;
        let cat_kwargs = vec![("sep", "\0")].into_py_dict(py);
        let strings: &PyString = str_accessor
            .call_method("cat", (), Some(cat_kwargs))?
            .downcast()?;
        Ok(Column {
            strings: strings
                .to_str()?
                .split('\0')
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        })
    }

    /// Iterator over the [Column]'s strings.
    pub fn iter(&self) -> impl Iterator<Item = &str> {
        self.strings.iter().copied()
    }
}

#[cfg(test)]
mod test_super {
    use super::*;
    use pyo3::types::{PyList, PyModule};
    use pyo3::Python;

    #[test]
    fn test_column() {
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
        let pandas = PyModule::import(py, "pandas").expect("can't import pandas");
        let arr = pandas
            .call1("Series", (PyList::new(py, vec!["foo", "bar", "baz"]),))
            .expect("can't construct a pandas Series");
        let column = Column::new(py, arr).expect("can't construct a Column");
        assert_eq!(column.iter().collect::<Vec<_>>(), vec!["foo", "bar", "baz"]);
    }
}
