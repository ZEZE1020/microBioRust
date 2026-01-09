**Core Infrastructure**
using Async where possible, parsing bioinformatic file types 
- [x] Genbank parser
- [x] Embl parser
- [x] Conversions of formats into gbk/embl/fasta/ffn/faa and save as gff/gbk or embl 
- [x] CI/CD with github actions and write tests - ongoing
- [x] Documentation with examples - ongoing

**Format expansion**
Integration of common types and parsers such as:
- [ ] VCF
- [x] BLAST output in various modes, -5 (XML) and -6 (one line per hit), compat with Diamond, MMSeqs2
- [ ] GTF
- [ ] link ups with formats from **rust-bio**, **noodles-vcf**
- [ ] look into metabolomics options (like KEGG) some are not open-source
- [ ] fastq gzipped version parsing
- [ ] SAM format parser
- [ ] RPKM output parser
- [ ] Further RNA-seq transcriptomics analysis
- [ ] Support for other compressed files such as BAM and CRAM
- [ ] Writer support for those types

**Advanced features**
Protein parameters such as: 
- [x] Hydrophobicity
- [x] Molecular weight
- [x]  Counted amino acids & as percentage
- [x] Aromaticity
- [ ]  Parsing phylogenetic trees 
- [x]  Parsing multiple sequence alignments
- [ ]  Methylation data parsing
- [x]  pyo3 integration

**Testing - bioinformatics**
- [ ] Testing & updating parsers with various edge case files such as gbk with the CONTIG(JOIN) structure

**Data Viz**
- [x] Example heatmap with WASM and Js example
- [ ]  Simple graph types bar, line, scatter
