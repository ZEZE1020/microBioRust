#[macro_export]
macro_rules! register_functions {
    ($module:expr, $($func:ident),*) => {
        $(
            // Using the full path pyo3::wrap_pyfunction ensures it always resolves
            $module.add_function(pyo3::wrap_pyfunction!($func, &$module)?)?;
        )*
    };
}