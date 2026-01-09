#![allow(unused_imports)]
use microBioRust::align::{load_msa_auto, AlignmentKind};
/// Tests for the multiple sequence alignment parsers
/// We can parse an MSA in fasta alignment format, iterate row-wise and column-wise
/// delete or subset the alignment and check consensus sequence or shannon entropy
/// The tests include a fasta format writer and print display in clustal format

// Unit tests
#[cfg(test)]
mod tests {
     use microBioRust::align::{load_msa_auto, AlignmentKind, Alignment, DNA, Protein};
     use std::path::PathBuf;

    // Helper to create a DNA alignment
    fn setup_test_dna_msa() -> Alignment<DNA> {
        let seq1 = b"ATGC--ATGCATGCATGC".to_vec(); // 3 rows of 6 columns
	let seq2 = b"ATGC-AATGCTTGCATGC".to_vec();
	let seq3 = b"TTGC-AATCCATGCAAGC".to_vec();
	let data : Vec<u8> = vec![seq1, seq2, seq3].into_iter().flatten().collect();
        let ids = vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()];
        Alignment::<DNA>::new(data, 3, ids)
    }
    #[test]
    #[allow(unused_mut)]
    fn test_msa_functions() {
        let msa = setup_test_dna_msa();
        let row_range = 0..msa.rows;
        let col_range = 0..msa.cols;
        let kind = AlignmentKind::DNA(msa);

        match kind {
            AlignmentKind::DNA(mut msa_inner) => {
                println!("Detected DNA alignment. Purging gaps > 50%...");
                msa_inner.purge_gappy_columns(0.5);
                let _ = msa_inner.display_interleaved(Some(row_range), Some(col_range), 60);
                msa_inner.write_fasta("cleaned_dna.fasta", None).unwrap();
            }
            AlignmentKind::Protein(mut msa_inner) => {
                println!("Detected Protein alignment. Calculating entropy...");
                let ent = msa_inner.column_entropy();
                let total: f32 = ent.into_iter().sum();
                println!("Column 0 Entropy: {:.2}", total);
                let _ = msa_inner.display_interleaved(Some(row_range), Some(col_range), 60);
                msa_inner.write_fasta("cleaned_protein.fasta", Some(60)).unwrap();
            }
         }
    }
}
