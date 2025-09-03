In microBioRust:

You can parse genbank files and save as a GFF (gff3) format as well as extracting DNA sequences, gene DNA sequences (ffn) and protein fasta sequences (faa)
Super simple way:

```rust
pub fn genbank_to_faa() -> Result<(), anyhow::Error> {
    let args = Arguments::parse();
    let records = genbank!(&args.filename);
    for record in records.iter() {
        for (k, v) in &record.cds.attributes {
            if let Some(seq) = record.seq_features.get_sequence_faa(k) {
                println!(">{}|{}\n{}", &record.id, &k, seq);
            }
        }
    }
    return Ok(());
}

```

Better for Debugging:

```rust
pub fn genbank_to_faa() -> Result<(), anyhow::Error> {
    let args: Vec<String> = env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem with parsing file arguments: {}", err);
        process::exit(1);
    });
    let file_gbk = fs::File::open(config.filename)?;
    let mut reader = Reader::new(file_gbk);
    let mut records = reader.records();
    let mut cds_counter: u32 = 0;
    loop {
        //collect from each record advancing on a next record basis, count cds records
        match records.next() {
            Some(Ok(mut record)) => {
                for (k, v) in &record.cds.attributes {
                    match record.seq_features.get_sequence_faa(&k) {
                        Some(value) => {
                            let seq_faa = value.to_string();
                            println!(">{}|{}\n{}", &record.id, &k, seq_faa);
                        }
                        _ => (),
                    };
                }
                cds_counter += 1;
            }
            Some(Err(e)) => {
                println!("Error encountered - an err {:?}", e);
            }
            None => {
                println!("finished iteration");
                break;
            }
        }
    }
    println!("Total records processed: {}", read_counter);
    return Ok(());
}
```

Example to save a provided multi- or single genbank file as a GFF file (by joining any multi-genbank)

```rust
pub fn genbank_to_gff() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let config = Config::new(&args).unwrap_or_else(|err| {
        println!("Problem with parsing file arguments: {}", err);
        process::exit(1);
        });
    let file_gbk = fs::File::open(&config.filename)?;
    let prev_start: u32 = 0;
    let mut prev_end: u32 = 0;
    let mut reader = Reader::new(file_gbk);
    let mut records = reader.records();
    let mut read_counter: u32 = 0;
    let mut seq_region: BTreeMap<String, (u32,u32)> = BTreeMap::new();
    let mut record_vec: Vec<Record> = Vec::new();
    loop {  
        match records.next() {   
            Some(Ok(mut record)) => {
                //println!("next record");
                //println!("Record id: {:?}", record.id);
                let source = record.source_map.source_name.clone().expect("issue collecting source name");
                let beginning = match record.source_map.get_start(&source) {
                                    Some(value) => value.get_value(),
                                    _ => 0,
                                    };
                let ending = match record.source_map.get_stop(&source) {
                                    Some(value) => value.get_value(),
                                    _ => 0,
                                    };
                if ending + prev_end < beginning + prev_end {
                    }
                seq_region.insert(source, (beginning + prev_end, ending + prev_end));
                record_vec.push(record);
                // Add additional fields to print if needed
                read_counter+=1;
                prev_end+=ending; // create the joined record if there are multiple
                },
            Some(Err(e)) => { println!("theres an err {:?}", e); },
                None => {
                println!("finished iteration");
                        break; },
            }
        }
    let output_file = format!("{}.gff", &config.filename);
    gff_write(seq_region.clone(), record_vec, &output_file, true);
    println!("Total records processed: {}", read_counter);
    return Ok(());
}
```

Example to create a completely new record, use of setters or set_ functionality

To write into GFF format requires gff_write(seq_region, record_vec, filename, true or false)

The seq_region is the region of interest to save with name and DNA coordinates such as `seqregion.entry("source_1".to_string(), (1,897))`

This makes it possible to save the whole file or to subset it

record_vec is a list of the records.  If there is only one record, include this as a vec using `vec![record]`

The boolean true/false describes whether the DNA sequence should be included in the GFF3 file

To write into genbank format requires gbk_write(seq_region, record_vec, filename), no true or false since genbank format will include the DNA sequence

 ```rust
pub fn create_new_record() -> Result<(), anyhow::Error> {
    let filename = format!("new_record.gff");
    let mut record = Record::new();
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    //example from E.coli K12
    seq_region.insert("source_1".to_string(), (1, 897));
    //Add the source into SourceAttributes
    record
        .source_map
        .set_counter("source_1".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(897))
        .set_organism("Escherichia coli".to_string())
        .set_mol_type("DNA".to_string())
        .set_strain("K-12 substr. MG1655".to_string())
        .set_type_material("type strain of Escherichia coli K12".to_string())
        .set_db_xref("PRJNA57779".to_string());
    //Add the features into FeatureAttributes, here we are setting two features, i.e. coding sequences or genes
    record
        .cds
        .set_counter("b3304".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(354))
        .set_gene("rplR".to_string())
        .set_product("50S ribosomal subunit protein L18".to_string())
        .set_codon_start(1)
        .set_strand(-1);
    record
        .cds
        .set_counter("b3305".to_string())
        .set_start(RangeValue::Exact(364))
        .set_stop(RangeValue::Exact(897))
        .set_gene("rplF".to_string())
        .set_product("50S ribosomal subunit protein L6".to_string())
        .set_codon_start(1)
        .set_strand(-1);
    //Add the sequences for the coding sequence (CDS) into SequenceAttributes
    record
        .seq_features
        .set_counter("b3304".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(354))
        .set_sequence_ffn(
            "ATGGATAAGAAATCTGCTCGTATCCGTCGTGCGACCCGCGCACGCCGCAAGCTCCAGGAG
CTGGGCGCAACTCGCCTGGTGGTACATCGTACCCCGCGTCACATTTACGCACAGGTAATT
GCACCGAACGGTTCTGAAGTTCTGGTAGCTGCTTCTACTGTAGAAAAAGCTATCGCTGAA
CAACTGAAGTACACCGGTAACAAAGACGCGGCTGCAGCTGTGGGTAAAGCTGTCGCTGAA
CGCGCTCTGGAAAAAGGCATCAAAGATGTATCCTTTGACCGTTCCGGGTTCCAATATCAT
GGTCGTGTCCAGGCACTGGCAGATGCTGCCCGTGAAGCTGGCCTTCAGTTCTAA"
                .to_string(),
        )
        .set_sequence_faa(
            "MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAE
QLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF"
                .to_string(),
        )
        .set_codon_start(1)
        .set_strand(-1);
    record
        .seq_features
        .set_counter("bb3305".to_string())
        .set_start(RangeValue::Exact(364))
        .set_stop(RangeValue::Exact(897))
        .set_sequence_ffn(
            "ATGTCTCGTGTTGCTAAAGCACCGGTCGTTGTTCCTGCCGGCGTTGACGTAAAAATCAAC
GGTCAGGTTATTACGATCAAAGGTAAAAACGGCGAGCTGACTCGTACTCTCAACGATGCT
GTTGAAGTTAAACATGCAGATAATACCCTGACCTTCGGTCCGCGTGATGGTTACGCAGAC
GGTTGGGCACAGGCTGGTACCGCGCGTGCCCTGCTGAACTCAATGGTTATCGGTGTTACC
GAAGGCTTCACTAAGAAGCTGCAGCTGGTTGGTGTAGGTTACCGTGCAGCGGTTAAAGGC
AATGTGATTAACCTGTCTCTGGGTTTCTCTCATCCTGTTGACCATCAGCTGCCTGCGGGT
ATCACTGCTGAATGTCCGACTCAGACTGAAATCGTGCTGAAAGGCGCTGATAAGCAGGTG
ATCGGCCAGGTTGCAGCGGATCTGCGCGCCTACCGTCGTCCTGAGCCTTATAAAGGCAAG
GGTGTTCGTTACGCCGACGAAGTCGTGCGTACCAAAGAGGCTAAGAAGAAGTAA"
                .to_string(),
        )
        .set_sequence_faa(
            "MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHADNTLTFGPRDGYAD
GWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSLGFSHPVDHQLPAG
ITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADEVVRTKEAKKK"
                .to_string(),
        )
        .set_codon_start(1)
        .set_strand(-1);
    //Add the full sequence of the entire record into the record.sequence
    record.sequence = "TTAGAACTGAAGGCCAGCTTCACGGGCAGCATCTGCCAGTGCCTGGACACGACCATGATA
TTGGAACCCGGAACGGTCAAAGGATACATCTTTGATGCCTTTTTCCAGAGCGCGTTCAGC
GACAGCTTTACCCACAGCTGCAGCCGCGTCTTTGTTACCGGTGTACTTCAGTTGTTCAGC
GATAGCTTTTTCTACAGTAGAAGCAGCTACCAGAACTTCAGAACCGTTCGGTGCAATTAC
CTGTGCGTAAATGTGACGCGGGGTACGATGTACCACCAGGCGAGTTGCGCCCAGCTCCTG
GAGCTTGCGGCGTGCGCGGGTCGCACGACGGATACGAGCAGATTTCTTATCCATAGTGTT
ACCTTACTTCTTCTTAGCCTCTTTGGTACGCACGACTTCGTCGGCGTAACGAACACCCTT
GCCTTTATAAGGCTCAGGACGACGGTAGGCGCGCAGATCCGCTGCAACCTGGCCGATCAC
CTGCTTATCAGCGCCTTTCAGCACGATTTCAGTCTGAGTCGGACATTCAGCAGTGATACC
CGCAGGCAGCTGATGGTCAACAGGATGAGAGAAACCCAGAGACAGGTTAATCACATTGCC
TTTAACCGCTGCACGGTAACCTACACCAACCAGCTGCAGCTTCTTAGTGAAGCCTTCGGT
AACACCGATAACCATTGAGTTCAGCAGGGCACGCGCGGTACCAGCCTGTGCCCAACCGTC
TGCGTAACCATCACGCGGACCGAAGGTCAGGGTATTATCTGCATGTTTAACTTCAACAGC
ATCGTTGAGAGTACGAGTCAGCTCGCCGTTTTTACCTTTGATCGTAATAACCTGACCGTT
GATTTTTACGTCAACGCCGGCAGGAACAACGACCGGTGCTTTAGCAACACGAGACAT"
        .to_string();
    gff_write(seq_region, vec![record], &filename, true);
    return Ok(());
}
```
