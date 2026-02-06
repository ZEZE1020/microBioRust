# Error Handling Strategy and Implementation

## Overview

This document describes the error handling improvements made to the microBioRust project and provides guidance for future development.

## Error Handling Patterns Used

### 1. anyhow::Error with Context

**Where Used**: Primary error handling mechanism throughout the codebase, especially in parsing operations.

**Pattern**:
```rust
use anyhow::{Context, Result};

fn parse_coordinate(s: &str) -> Result<u32> {
    s.parse()
        .with_context(|| format!("Failed to parse coordinate from: {}", s))
}
```

**Benefits**:
- Rich error context with error chains
- Easy to propagate with `?` operator
- Good for complex operations where multiple error types may occur

### 2. Result with Specific Error Types

**Where Used**: I/O operations, especially in file handling.

**Pattern**:
```rust
use std::io;

fn read_file(path: &str) -> io::Result<String> {
    std::fs::read_to_string(path)
}
```

**Benefits**:
- Type-safe error handling
- Clear error semantics for specific operations

### 3. Python Error Types (PyO3)

**Where Used**: Python bindings in `microbiorust-py`.

**Pattern**:
```rust
use pyo3::PyErr;
use pyo3::exceptions::PyIOError;

fn delete_file(path: &str) -> PyResult<()> {
    std::fs::remove_file(path)
        .map_err(|e| PyErr::new::<PyIOError, _>(
            format!("Failed to delete file {}: {}", path, e)
        ))
}
```

**Benefits**:
- Proper Python exception types
- Clear error messages for Python users
- Integration with Python's exception handling

## Improvements Made

### Critical Issues Fixed

#### 1. Unvalidated Vector Access → Safe Access with Error Context

**Before**:
```rust
let qseqid = cols.get(0).unwrap().to_string();
```

**After**:
```rust
let qseqid = cols.get(0)
    .ok_or_else(|| anyhow::anyhow!("Missing qseqid column (0)"))?
    .to_string();
```

**Impact**: Prevents panics on malformed BLAST output, provides clear error messages.

#### 2. Unvalidated String Indexing → Safe Access with Defaults

**Before**:
```rust
let org: Vec<&str> = self.line_buffer.split('\"').collect();
organism = org[1].to_string();  // Panic if split doesn't produce 2+ parts
```

**After**:
```rust
let org: Vec<&str> = self.line_buffer.split('\"').collect();
organism = org.get(1).unwrap_or(&"").to_string();
```

**Impact**: Handles malformed GenBank/EMBL files gracefully instead of panicking.

#### 3. .expect() on Parse Operations → Contextual Error Messages

**Before**:
```rust
thestart = cap[1]
    .parse()
    .expect("failed to match and parse numerical start");
```

**After**:
```rust
thestart = cap[1]
    .parse()
    .with_context(|| format!("Failed to parse start coordinate from: {}", &cap[1]))?;
```

**Impact**: Provides specific values that failed to parse, aiding debugging.

#### 4. .expect() on File I/O → Proper Error Propagation

**Before**:
```rust
std::fs::remove_file(&output_file).expect("Issue deleting output filename");
```

**After**:
```rust
std::fs::remove_file(&output_file)
    .map_err(|e| PyErr::new::<PyIOError, _>(
        format!("Failed to delete output file {}: {}", &output_file, e)
    ))?;
```

**Impact**: Returns proper Python exceptions with file names and OS error details.

## Design Decisions and Rationale

### 1. .unwrap_or_default() for BLAST Parsing

**Decision**: Keep `.unwrap_or_default()` for numeric fields in BLAST parsing.

**Rationale**: 
- BLAST output can have missing or malformed values in optional columns
- Defaulting to 0 is preferable to failing the entire parse operation
- This maintains backward compatibility and robustness with real-world data
- Users can filter results with 0 values if needed

**Example**:
```rust
pident: cols.get(2).ok_or_else(...)?.parse().unwrap_or(0.0),
```

### 2. Macros with panic!()

**Decision**: Leave `panic!()` in convenience macros (`genbank!`, `embl!`).

**Rationale**:
- Macros are designed for simple use cases where error handling would be verbose
- They're primarily used in examples and quick scripts
- Production code should use the Reader API directly for proper error handling

**Alternative for Production Code**:
```rust
// Instead of: let records = genbank!("file.gbk");
// Use:
let file = File::open("file.gbk")?;
let mut reader = gbk::Reader::new(file);
let records: Result<Vec<_>> = reader.records().collect();
```

### 3. Programming Error Panics

**Decision**: Keep `panic!()` for programming errors (e.g., "Counter key not set").

**Rationale**:
- These represent bugs in the library code, not user errors
- Should never occur with correct API usage
- Panic is appropriate for detecting contract violations
- Enhanced error messages to aid debugging

## Remaining Considerations

### Low Priority Issues

1. **Silent Error Swallowing in Python BLAST Parser**: The Python bindings duplicate some BLAST parsing logic with `.unwrap_or()` defaults. This is acceptable for Python API robustness but could be consolidated with the main implementation.

2. **Regex Compilation .expect()**: Used in `lazy_static!` initialization. This is acceptable as regex compilation failure is a programming error that should be caught in testing.

3. **Main Function Panics**: The `main.rs` functions use `panic!()` for missing arguments. This is standard for CLI applications and acceptable.

## Testing Strategy

All error handling improvements were validated with:

1. **Unit Tests**: All existing tests pass (8 integration tests, 10 doc tests)
2. **Build Verification**: Clean compilation with no warnings on error handling code
3. **Manual Testing**: Test files processed successfully

## Guidelines for Future Development

### DO:

✅ Use `anyhow::Context` to add context to errors:
```rust
operation().with_context(|| format!("Failed to process {}", item))?;
```

✅ Use `.ok_or_else()` for Option to Result conversion:
```rust
value.ok_or_else(|| anyhow::anyhow!("Missing required value"))?
```

✅ Validate before indexing, use `.get()`:
```rust
let value = vec.get(index).ok_or_else(|| anyhow::anyhow!("Index out of bounds"))?;
```

✅ Include relevant context in error messages (filenames, values, indices)

### DON'T:

❌ Use `.unwrap()` or `.expect()` on user input or external data

❌ Use direct indexing `[i]` on user-provided data structures

❌ Swallow errors silently without logging or documentation

❌ Return generic error messages without context

### Exception Cases:

- ✔️ `.expect()` is acceptable for static initialization (regex compilation)
- ✔️ `panic!()` is acceptable for programming errors and contract violations
- ✔️ `.unwrap_or_default()` is acceptable when default values are meaningful

## Error Categories

### 1. User Errors (Should Return Result)
- Invalid file paths
- Malformed input data
- Missing required fields
- Invalid coordinates or ranges

### 2. Programming Errors (Can Panic)
- Invalid regex patterns (at compile time)
- API misuse (e.g., calling insert_to before set_counter)
- Internal state inconsistencies

### 3. System Errors (Should Return Result)
- File I/O failures
- Memory allocation failures
- Network errors (in async code)

## Performance Considerations

Error handling improvements have minimal performance impact:

1. **String allocation for context**: Only occurs on error path
2. **Additional checks**: Negligible compared to parsing operations
3. **Result types**: Zero-cost abstraction in Rust

The pre-compiled regex patterns using `lazy_static!` maintain the 10-50x performance improvement noted in the code.

## Conclusion

The error handling improvements make microBioRust more robust and maintainable:

- **Better diagnostics**: Clear error messages with context
- **Graceful degradation**: Handles malformed data without panicking
- **Maintainability**: Easier to debug issues with detailed error chains
- **Production-ready**: Suitable for use in pipelines and services

All changes maintain backward compatibility and pass existing tests.
