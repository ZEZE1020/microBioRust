import pytest
import os
import microbiorust
import textwrap
from pathlib import Path

# --- FIXTURES: Generating Mock Data ---

@pytest.fixture
def mock_gbk(tmp_path):
    """Creates a minimal valid GenBank file."""
    path = tmp_path / "test.gbk"
    content = textwrap.dedent("""
        LOCUS       source_1                 910 bp    DNA     linear   CON 01-NOV-2024
        DEFINITION  Escherichia coli K-12 substr. MG1655.
        ACCESSION   source_1
        VERSION     source_1
        KEYWORDS    .
        SOURCE      Escherichia coli K-12 substr. MG1655
          ORGANISM  Escherichia coli K-12 substr. MG1655
        FEATURES             Location/Qualifiers
             source          1..910
                             /organism="K-12 substr. MG1655"
                             /mol_type="DNA"
             gene            complement(1..354)
                             /locus_tag="b3304"
             CDS             complement(1..354)
                             /locus_tag="b3304"
                             /codon_start=1
                             /gene="rplR"
                             /translation="MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALDAAREAGLQ"
                             /product="50S ribosomal subunit protein L18"
             gene            complement(364..897)
                             /locus_tag="b3305"
             CDS             complement(364..897)
                             /locus_tag="b3305"
                             /codon_start=1
                             /gene="rplF"
                             /translation="MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHNTLTFGPRDGYADGWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSGFSHPVDHQLPAGITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADVVRTKEAKK"
                             /product="50S ribosomal subunit protein L6"
        ORIGIN
                1 TTAGAACTGA AGGCCAGCTT CACGGGCAGC ATCTGCCAGT GCCTGGACAC GACCATGATA
               61 TTGGAACCCG GAACGGTCAA AGGATACATC TTTGATGCCT TTTTCCAGAG CGCGTTCAGC
              121 GACAGCTTTA CCCACAGCTG CAGCCGCGTC TTTGTTACCG GTGTACTTCA GTTGTTCAGC
              181 GATAGCTTTT TCTACAGTAG AAGCAGCTAC CAGAACTTCA GAACCGTTCG GTGCAATTAC
              241 CTGTGCGTAA ATGTGACGCG GGGTACGATG TACCACCAGG CGAGTTGCGC CCAGCTCCTG
              301 GAGCTTGCGG CGTGCGCGGG TCGCACGACG GATACGAGCA GATTTCTTAT CCATAGTGTT
              361 ACCTTACTTC TTCTTAGCCT CTTTGGTACG CACGACTTCG TCGGCGTAAC GAACACCCTT
              421 GCCTTTATAA GGCTCAGGAC GACGGTAGGC GCGCAGATCC GCTGCAACCT GGCCGATCAC
              481 CTGCTTATCA GCGCCTTTCA GCACGATTTC AGTCTGAGTC GGACATTCAG CAGTGATACC
              541 CGCAGGCAGC TGATGGTCAA CAGGATGAGA GAAACCCAGA GACAGGTTAA TCACATTGCC
              601 TTTAACCGCT GCACGGTAAC CTACACCAAC CAGCTGCAGC TTCTTAGTGA AGCCTTCGGT
              661 AACACCGATA ACCATTGAGT TCAGCAGGGC ACGCGCGGTA CCAGCCTGTG CCCAACCGTC
              721 TGCGTAACCA TCACGCGGAC CGAAGGTCAG GGTATTATCT GCATGTTTAA CTTCAACAGC
              781 ATCGTTGAGA GTACGAGTCA GCTCGCCGTT TTTACCTTTG ATCGTAATAA CCTGACCGTT
              841 GATTTTTACG TCAACGCCGG CAGGAACAAC GACCGGTGCT TTAGCAACAC GAGACA
        //
    """)
    path.write_text(content)
    return str(path)

@pytest.fixture
def mock_msa(tmp_path):
    """Creates a mock FASTA alignment file."""
    path = tmp_path / "align.fasta"
    content = (
        ">Seq1\nATGC--AT\n"
        ">Seq2\nATGC--TT\n"
        ">Seq3\nATGCGGTT\n"
    )
    path.write_text(content)
    return str(path)

@pytest.fixture
def mock_blast_tab(tmp_path):
    """Creates a mock BLAST tabular file."""
    path = tmp_path / "results.tab"
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    content = "seqA\tseqB\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-10\t200.0\n"
    path.write_text(content)
    return str(path)

#tests for parsers (gbk, embl)

def test_gbk_conversions(mock_gbk):
    #test protein sequence extraction
    #with open(mock_gbk, 'rb') as f:
    #    print(f"\n[DEBUG] Raw File Content: {f.read()}")
    faa = microbiorust.gbk_to_faa(mock_gbk)
    assert os.path.exists(mock_gbk)
    assert any("b3304" in line for line in faa)
    
    #test count
    count = microbiorust.gbk_to_faa_count(mock_gbk)
    assert count == 2
    
    #test nucleotide sequence extraction
    fna = microbiorust.gbk_to_fna(mock_gbk)
    assert any("ttagaactga" in line.lower() for line in fna)

def test_gbk_to_gff(mock_gbk):
    #this function writes to {filename}.gff and reads again to check
    microbiorust.gbk_to_gff(mock_gbk, dna=True)
    gff_path = f"{mock_gbk}.gff"
    assert os.path.exists(gff_path)
    with open(gff_path, "r") as f:
        assert "source_source_1_1" in f.read()

#tests for multiple sequence alignment 

def test_subset_msa(mock_msa):
    #subset mock alignment: Rows 0-2 (Seq1 & Seq2), Cols 0-4 (ATGC)
    subset = microbiorust.subset_msa_alignment(mock_msa, (0, 2), (0, 4))
    print("subset", subset)
    assert len(subset) == 2  # 2 headers + 2 sequences
    assert ">Seq1" in subset[0]
    assert "ATGC" in subset[1]

def test_purge_gaps(mock_msa, tmp_path):
    out_path = str(tmp_path / "purged.fasta")
    #threshold 0.5 should remove the '--' columns in Seq1 and Seq2 and write to file
    microbiorust.purge_gaps(mock_msa, out_path, 0.5)
    assert os.path.exists(out_path)

def test_get_consensus(mock_msa):
    #given 'ATGC' is constant in the mock, it should be in consensus
    consensus = microbiorust.get_consensus(mock_msa)
    assert consensus.startswith("ATGC")

#tests for Sequence Metrics

def test_hydrophobicity():
    seq = "MALWMRLLPLLALLALWGPDPAAAFVN"
    scores = microbiorust.hydrophobicity(seq, window_size=3)
    assert len(scores) > 0
    assert all(isinstance(s, float) for s in scores)

def test_amino_counts():
    seq = "MATAG"
    counts = microbiorust.amino_counts(seq)
    assert counts['M'] == 1
    assert counts['A'] == 2

#test for Async Tabular Parser

def test_parse_tabular(mock_blast_tab):
    results = microbiorust.parse_tabular(mock_blast_tab)
    assert len(results) == 1
    assert results[0]['qseqid'] == "seqA"
    assert results[0]['bitscore'] == 200.0
