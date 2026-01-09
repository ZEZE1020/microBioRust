//! A Multiple Sequence Alignment parser (MSA)
//!
//! You are able to parse multiple sequence alignments in either DNA or protein as fasta alignment format
//! You can collect a row by index or by ID
//! It is also possible to delete a row by an indexed number (and subsequently reindex)
//! You can count columns, identify columns with gaps, or subset the alignment
//! It is also possible to iterate row-wise or column-wise, identify base frequencies at a particular column
//! measure shannon entropy or identify a consensus.  You can also print either the full or partial alignment
//!
//!
//! use microBioRust::align::{load_msa_auto, AlignmentKind};
//! use clap::Parser;
//! use std::path::PathBuf;
//!
//! #[derive(Parser, Debug)]
//! #[command(author, version, about = "microBioRust: Sustainable MSA Processing Tool")]
//! struct Args {
//!    /// Path to the input FASTA/MSA file
//!    input: PathBuf,
//!    ///number of rows to display (e.g., --rows 10)
//!    #[arg(short, long)]
//!    rows: Option<usize>,
//!    //number of columns to display (e.g., --cols 100)
//!    #[arg(short, long)]
//!    cols: Option<usize>,
//!    //width of the interleaved blocks (default: 60)
//!    #[arg(short, long, default_value_t = 60)]
//!    width: usize,
//! }
//!
//! #[allow(unused_mut)]
//! fn main() {
//!    let args = Args::parse();
//!    println!("-- microBioRust v{} ---", env!("CARGO_PKG_VERSION"));
//!
//!
//! let msa_kind = load_msa_auto(&args.input).expect("file not found error");
//!    let row_range = args.rows.map(|r| 0..r);
//!    let col_range = args.cols.map(|c| 0..c);
//!    
//!   match Ok::<AlignmentKind, anyhow::Error>(msa_kind) {
//!        Ok(AlignmentKind::DNA(mut msa)) => {
//!            println!("Detected DNA alignment. Purging gaps > 50%...");
//!            msa.purge_gappy_columns(0.5);
//!            msa.display_interleaved(row_range, col_range, args.width);
//!            msa.write_fasta("cleaned_dna.fasta", None).unwrap();
//!        }
//!        Ok(AlignmentKind::Protein(mut msa)) => {
//!            println!("Detected Protein alignment. Calculating entropy...");
//!            let ent = msa.column_entropy();
//!            let total: f32 = ent.into_iter().sum();
//!            println!("Column 0 Entropy: {:.2}", total);
//!            msa.display_interleaved(row_range, col_range, args.width);
//!            msa.write_fasta("cleaned_protein.fasta", Some(60)).unwrap();
//!        }
//!        Err(e) => eprintln!("Error: {}", e),
//!    }
//! }

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write, Result};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::ops::Range;
use anyhow::{anyhow, Context};


pub struct DNA;
pub struct Protein;

pub struct Alignment<Kind> {
    data: Vec<u8>,
    num_rows: usize,
    num_cols: usize,
    ids: Vec<String>,
    id_to_index: HashMap<String, usize>,
    _marker: std::marker::PhantomData<Kind>,
}

impl<Kind> Alignment<Kind> {
    pub fn new(data: Vec<u8>, num_rows: usize, ids: Vec<String>) -> Self {
        let num_cols = data.len() / num_rows;
        let id_to_index = ids.iter().enumerate().map(|(i, id)| (id.clone(), i)).collect();
        Self { data, num_rows, num_cols, ids, id_to_index, _marker: std::marker::PhantomData }
    }
    //to access a row by ID
    pub fn get_row_by_id(&self, id: &str) -> Option<&[u8]> {
        let &idx = self.id_to_index.get(id)?;
        let start = idx * self.num_cols;
        Some(&self.data[start..start + self.num_cols])
    }
    //to delete a column by start and end
    pub fn remove_columns(&mut self, start: usize, end: usize) {
        let cols_to_remove = end - start;
        let mut new_data = Vec::with_capacity(self.data.len() - (self.num_rows * cols_to_remove));

        for row in self.data.chunks_exact(self.num_cols) {
            //keeping everything before the slice and everything after
            new_data.extend_from_slice(&row[..start]);
            new_data.extend_from_slice(&row[end..]);
        }

        self.data = new_data;
        self.num_cols -= cols_to_remove;
    }
    //to delete a row by indexed number
    pub fn remove_row_by_number(&mut self, row_idx: usize) {
        if row_idx >= self.num_rows { return; }
        
        let start = row_idx * self.num_cols;
        let end = start + self.num_cols;
        
        // Remove from the data vector
        self.data.drain(start..end);
        
        // Update the ID map and row count
        let id_to_remove = self.ids[row_idx].clone();
        self.id_to_index.remove(&id_to_remove);
        self.ids.remove(row_idx);
        self.num_rows -= 1;
        
        // Re-index the map (crucial for maintaining O(1) lookup)
        self.reindex_map();
    }
    //to reindex the map (if items have been deleted
    fn reindex_map(&mut self) {
        for (i, id) in self.ids.iter().enumerate() {
            self.id_to_index.insert(id.clone(), i);
        }
    }
    //to remove a row by an ID string as &str
    pub fn remove_row_by_id(&mut self, id: &str) {
        if let Some(&idx) = self.id_to_index.get(id) {
            self.remove_row_by_number(idx);
            }
    }
    //to get a frequency map of residues at a given column index.
    pub fn column_counts(&self, col_idx: usize) -> HashMap<u8, usize> {
        let mut counts = HashMap::new();
        
        for r in 0..self.num_rows {
            let residue = self.data[r * self.num_cols + col_idx];
            // The .entry() API is the "Senior" way to do this in Rust
            *counts.entry(residue).or_insert(0) += 1;
        }
        
        counts
    }
    //to get a list vector of counters for every column in the alignment
    pub fn all_column_counts(&self) -> Vec<HashMap<u8, usize>> {
        (0..self.num_cols)
            .map(|c| self.column_counts(c))
            .collect()
    }
    //to check for gaps - Returns true if the proportion of gaps '-' in a column exceeds the provided cutoff (0.0 to 1.0)
    pub fn is_gap_heavy(&self, col_idx: usize, cutoff: f32) -> bool {
        let mut gap_count = 0;
        
        for r in 0..self.num_rows {
            if self.data[r * self.num_cols + col_idx] == b'-' {
                gap_count += 1;
            }
        }
        
        let gap_proportion = gap_count as f32 / self.num_rows as f32;
        gap_proportion > cutoff
    }
    //to identify columns that exceed the provided gap threshold
    //returns a list of indices to be removed
    pub fn identify_gappy_columns(&self, cutoff: f32) -> Vec<usize> {
        (0..self.num_cols)
            .filter(|&c| self.is_gap_heavy(c, cutoff))
            .collect()
    }
    pub fn purge_gappy_columns(&mut self, cutoff: f32) {
        let to_remove = self.identify_gappy_columns(cutoff);
        
        //NB: when removing multiple columns, we must remove from right to left to avoid index shifting, or use a new buffer
        for &col_idx in to_remove.iter().rev() {
            self.remove_columns(col_idx, col_idx + 1);
        }
    }
    //to print the alignment. 
    //row_range: e.g., Some(0..5) for first five rows. None for all.
    //col_range: e.g., Some(0..50) for first 50 columns. None for all.
    pub fn display(&self, row_range: Option<Range<usize>>, col_range: Option<Range<usize>>) -> io::Result<()> {
        // Fallback to full range if None is provided
        let rows = row_range.unwrap_or(0..self.num_rows);
        let cols = col_range.unwrap_or(0..self.num_cols);

        let stdout = io::stdout();
        let mut handle = BufWriter::new(stdout.lock());

        // Iterate only through the requested rows
        for r_idx in rows {
            // Safety check: ensure the index is within bounds
            if let Some(id) = self.ids.get(r_idx) {
                writeln!(handle, ">{}", id)?;

                // MATH: Find the start of the row, then offset by the column start
                let row_start_in_buffer = r_idx * self.num_cols;
                let start = row_start_in_buffer + cols.start;
                let end = row_start_in_buffer + cols.end;

                // Ensure we don't slice past the end of the actual row
                let safe_end = end.min(row_start_in_buffer + self.num_cols);
                
                if let Some(row_slice) = self.data.get(start..safe_end) {
                    handle.write_all(row_slice)?;
                    handle.write_all(b"\n")?;
                }
            }
        }
        
        handle.flush()?;
        Ok(())
    }
    pub fn display_interleaved(
        &self, 
        row_range: Option<Range<usize>>, 
        col_range: Option<Range<usize>>,
        block_width: usize // Usually 60 or 80
    ) -> io::Result<()> {
        let rows = row_range.unwrap_or(0..self.num_rows);
        let cols = col_range.unwrap_or(0..self.num_cols);

        let stdout = io::stdout();
        let mut handle = BufWriter::new(stdout.lock());
        let max_id_width = rows.clone().map(|i| self.ids[i].len()).max().unwrap_or(10).min(20); //capped at 20
        // 1. Iterate through columns in "blocks"
        for block_start in (cols.start..cols.end).step_by(block_width) {
            let block_end = (block_start + block_width).min(cols.end);

            // 2. For each block, print every requested row
            for r_idx in rows.clone() {
                if let Some(id) = self.ids.get(r_idx) {
                    //truncate or pad ID for alignment
                    let display_id = if id.len() > 10 { &id[..10] } else { id };
                    
                    // print ID with padding for clean columns
                    write!(handle, "{:<width$} ", display_id, width=max_id_width)?;

                    //slice the specific chunk for this row
                    let row_offset = r_idx * self.num_cols;
                    let start = row_offset + block_start;
                    let end = row_offset + block_end;

                    if let Some(slice) = self.data.get(start..end) {
                        handle.write_all(slice)?;
                        writeln!(handle)?;
                    }
                }
            }
            //add a spacer between blocks
            writeln!(handle)?;
        }
        
        handle.flush()?;
        Ok(())
    }

    //to subset the alignment if that is needed - by row and column
    pub fn subset(&self, rows: Range<usize>, cols: Range<usize>) -> Self {
        let new_num_rows = rows.len();
        let new_num_cols = cols.len();
        let mut new_data = Vec::with_capacity(new_num_rows * new_num_cols);
        let mut new_ids = Vec::with_capacity(new_num_rows);

        for r in rows {
            // 1. Grab the ID
            new_ids.push(self.ids[r].clone());

            // 2. Grab the specific column slice for this row
            let start = (r * self.num_cols) + cols.start;
            let end = (r * self.num_cols) + cols.end;
            new_data.extend_from_slice(&self.data[start..end]);
        }

        Self::new(new_data, new_num_rows, new_ids)
    }
    //to write part of the alignment to file
    pub fn write_fasta_part<P: AsRef<Path>>(
        &self, 
        path: P, 
        row_range: Option<Range<usize>>, 
        col_range: Option<Range<usize>>,
        wrap: Option<usize>
    ) -> anyhow::Result<()> {
        let rows = row_range.unwrap_or(0..self.num_rows);
        let cols = col_range.unwrap_or(0..self.num_cols);
        
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        for r in rows {
            writeln!(writer, ">{}", self.ids[r])?;
            
            let start = (r * self.num_cols) + cols.start;
            let end = (r * self.num_cols) + cols.end;
            let row_chunk = &self.data[start..end];

            if let Some(w) = wrap {
                for chunk in row_chunk.chunks(w) {
                    writer.write_all(chunk)?;
                    writer.write_all(b"\n")?;
                }
            } else {
                writer.write_all(row_chunk)?;
                writer.write_all(b"\n")?;
            }
        }
        Ok(())
    }
    //to write the current state of the alignment to a FASTA file.
    //to wrap: Optional number of characters per line (e.g., Some(60)). 
    //or use None for a single line per sequence.
    pub fn write_fasta<P: AsRef<Path>>(&self, path: P, wrap: Option<usize>) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        for r in 0..self.num_rows {
            // Write the Header line
            writeln!(writer, ">{}", self.ids[r])?;

            let start = r * self.num_cols;
            let end = start + self.num_cols;
            let row_data = &self.data[start..end];

            match wrap {
                Some(width) => {
                    // Write in wrapped chunks for legacy compatibility
                    for chunk in row_data.chunks(width) {
                        writer.write_all(chunk)?;
                        writer.write_all(b"\n")?;
                    }
                }
                None => {
                    // Modern single-line writing
                    writer.write_all(row_data)?;
                    writer.write_all(b"\n")?;
                }
            }
        }
        
        //ensure all the data is flushed to the disk
        writer.flush()?;
        Ok(())
    }
}
impl Alignment<DNA> {
    //calculate GC-content of a specific column (homologous site)
    pub fn gc_content_at_col(&self, col_idx: usize) -> f32 {
        let mut gc_count = 0;
        for r in 0..self.num_rows {
            let base = self.data[r * self.num_cols + col_idx].to_ascii_uppercase();
            if base == b'G' || base == b'C' {
                gc_count += 1;
            }
        }
        gc_count as f32 / self.num_rows as f32
    }
    //to determine the most frequent base at a column
    pub fn most_frequent_base(&self, column: &[u8]) -> char {
        let mut counts = HashMap::new();
        for &base in column.iter().filter(|&&b| b != b'-') {
            *counts.entry(base).or_insert(0) += 1;
        }
        
        counts.into_iter()
            .max_by_key(|&(_, count)| count)
            .map(|(base, _)| base as char)
            .unwrap_or('-') // Default to gap if column is empty
    }
    //to provide a consensus sequence
    pub fn consensus_sequence(&self) -> String {
        self.iter_cols()
            .map(|col| {
                // Find the most frequent base in this column
                // Senior move: ignore gaps '-' during consensus calculation
                self.most_frequent_base(&col)
            })
            .collect()
    }
    //to identify a degenerate site
    pub fn is_degenerate_site(&self, col_idx: usize, threshold: f32) -> bool {
        let counts = self.column_counts(col_idx);
        let max_freq = *counts.values().max().unwrap_or(&0) as f32 / self.num_rows as f32;
        max_freq < threshold
    }
}
//for protein MSA
impl Alignment<Protein> {
    //to find a single shannon entropy value for a column
    pub fn calculate_shannon_entropy(&self, column: &[u8]) -> f32 {
        if column.is_empty() { return 0.0; }
        let mut counts = HashMap::new();
	let mut non_gap_count = 0;
	//only counting non gaps
	for &residue in column.iter().filter(|&&b| b != b'-') {
	  *counts.entry(residue).or_insert(0) += 1;
	  non_gap_count +=1;
	  }
	if non_gap_count == 0 { return 0.0; }
	let total = non_gap_count as f32;
        counts.values()
            .map(|&count| {
                let p = count as f32/total;
		-p * p.log2()
		}).sum()
            }
    //to find the column entropy of a group of columns as a list
    pub fn column_entropy(&self) -> Vec<f32> {
        self.iter_cols().map(|col| self.calculate_shannon_entropy(&col)).collect()
	}
    //to identify rare residues at a column index
    pub fn rare_residues_at_site(&self, col_idx: usize) -> Vec<u8> {
        let counts = self.column_counts(col_idx);
        counts.into_iter()
            .filter(|&(_, count)| count == 1) // Residues appearing only once
            .map(|(residue, _)| residue)
            .collect()
    }
}
//to iterate column-wise
pub struct ColumnIterator<'a, Kind> {
    alignment: &'a Alignment<Kind>,
    current_col: usize,
}

impl<'a, Kind> Iterator for ColumnIterator<'a, Kind> {
    type Item = Vec<u8>; // The residues at this homologous site

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_col >= self.alignment.num_cols {
            return None;
        }

        // Collect residues from each row for the current column
        let column: Vec<u8> = (0..self.alignment.num_rows)
            .map(|r| {
                let idx = (r * self.alignment.num_cols) + self.current_col;
                self.alignment.data[idx]
            })
            .collect();

        self.current_col += 1;
        Some(column)
    }
}

//helper to the main Alignment struct
impl<Kind> Alignment<Kind> {
    pub fn iter_cols(&self) -> ColumnIterator<'_, Kind> {
        ColumnIterator {
            alignment: self,
            current_col: 0,
        }
    }
}
//enum to determine if the MSA is DNA or protein
pub enum AlignmentKind {
    DNA(Alignment<DNA>),
    Protein(Alignment<Protein>),
}
//to load the MSA we first parse the fasta from path
pub fn parse_fasta<P: AsRef<Path>>(filename: P) -> anyhow::Result<(Vec<String>, Vec<u8>)> {
    let file = File::open(&filename).context("Failed to open FASTA file")?;
    let reader = BufReader::new(file);
    let mut found_first_header = false;
    let mut ids = Vec::new();
    let mut sequences = Vec::new();
    let mut current_seq = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read line")?;
        let line = line.trim();
        
        if line.is_empty() || line.starts_with('#') { continue; }

        if line.starts_with('>') {
            found_first_header = true;
            if !current_seq.is_empty() {
                sequences.push(current_seq);
                current_seq = Vec::new();
            }
            // Store the ID (minus the '>')
            ids.push(line[1..].to_string());
        } else {
	    if !found_first_header {
	       return Err(anyhow!("format error line {}: please ensure the MSA is in fasta format (headers with >)", &line_idx));
	       }
            //append sequence data (handling multi-line sequences)
            current_seq.extend_from_slice(line.as_bytes());
        }
    }
    // Push the final sequence
    if !current_seq.is_empty() {
        sequences.push(current_seq);
    }

    //to validate: Ensure all sequences are the same length (required for MSA)
    let first_len = sequences.first().map(|s| s.len()).unwrap_or(0);
    if !sequences.iter().all(|s| s.len() == first_len) {
        return Err(anyhow!("validation error: Sequences have unequal lengths."));
    }

    // Flatten into a single Vec<u8> for our Alignment struct
    let flat_data = sequences.into_iter().flatten().collect();

    Ok((ids, flat_data))
}
//to load the MSA
pub fn load_msa_auto<P: AsRef<Path>>(filename: P) -> anyhow::Result<AlignmentKind> {
    //read file
    let (ids, data) = parse_fasta(&filename)?;
    if ids.is_empty() || data.is_empty() {
        return Err(anyhow!("empty input cannot determine alignment type"));
	}
    //automagic detection
    let is_protein = data.iter()
        .filter(|&&b| b != b'-') // Ignore gaps
        .take(100)               // Peek at first 100 residues
        .any(|&b| {
            let c = b.to_ascii_uppercase();
            // Amino acids that don't exist in DNA/RNA
            matches!(c, b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' | b'R' | b'V' | b'W' | b'Y')
        });

    let num_rows = ids.len();
    if is_protein {
        Ok(AlignmentKind::Protein(Alignment::<Protein>::new(data, num_rows, ids)))
    } else {
        Ok(AlignmentKind::DNA(Alignment::<DNA>::new(data, num_rows, ids)))
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a dummy DNA alignment
    fn setup_test_dna_msa() -> Alignment<DNA> {
        let seq1 = b"ATGC--ATGCATGCATGC".to_vec(); // 3 rows of 6 columns
	let seq2 = b"ATGC-AATGCTTGCATGC".to_vec();
	let seq3 = b"TTGC-AATCCATGCAAGC".to_vec();
	let data : Vec<u8> = vec![seq1, seq2, seq3].into_iter().flatten().collect();
        let ids = vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()];
        Alignment::<DNA>::new(data, 3, ids)
    }

    #[test]
    fn test_row_lookup_by_id() {
        let msa = setup_test_dna_msa();
        let row = msa.get_row_by_id("seq2").expect("ID should exist");
        assert_eq!(row, b"ATGC-AATGCTTGCATGC"); // Check the middle row
    }

    #[test]
    fn test_column_iteration_and_counts() {
        let msa = setup_test_dna_msa();
        //column 0 is A, A, T
        let counts = msa.column_counts(0);
        assert_eq!(*counts.get(&b'A').unwrap(), 2);
        assert_eq!(*counts.get(&b'T').unwrap(), 1);
    }

    #[test]
    fn test_gap_heavy_detection() {
        let msa = setup_test_dna_msa();
        //in this example column 5 has 1 gap out of 3 (33%)
        assert!(msa.is_gap_heavy(5, 0.3));   // 0.33 > 0.3 is true
        assert!(!msa.is_gap_heavy(5, 0.4));  // 0.33 > 0.4 is false
    }

    #[test]
    fn test_column_removal() {
        let mut msa = setup_test_dna_msa();
        let original_cols = msa.num_cols;
    
        // Remove columns 1 and 2 (indices 1..3)
        msa.remove_columns(1, 3);
    
        assert_eq!(msa.num_cols, original_cols - 2);
    
        // row 0 was: A [T G] C - - A T G C A T G C A T G C
        // should be: A C - - A T G C A T G C A T G C
        let result = msa.get_row_by_id("seq1").unwrap();
        assert_eq!(result, b"AC--ATGCATGCATGC");
   }

   #[test]
   fn test_protein_column_entropy() {
       // Column 0: Perfectly conserved (M, M) -> Entropy 0
       // Column 1: Half M, half A -> Entropy 1.0
       let data = b"MMMA".to_vec(); 
       let ids = vec!["p1".into(), "p2".into()];
       let msa = Alignment::<Protein>::new(data, 2, ids);
       let result = msa.column_entropy();
       let total: f32 = result.iter().sum();
       assert!((total - 1.0).abs() < f32::EPSILON);
       }
}
