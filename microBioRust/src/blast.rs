#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(unused_variables)]
#![allow(dead_code)]
//! async streaming BLAST parsers for the (outfmt 6 and outfmt 5),
//! gzip-transparent reader, option of --json streaming output
//!
//! Async streaming BLAST parsers
//! - outfmt 6 (tabular): streamed, async line reader -> JSON NDJSON
//! - outfmt 5 (XML): async streaming XML iterator returning Iteration structs
//! - Gzip auto-detection by filename (".gz")
//! - CLI with --json to stream newline-delimited JSON to stdout
//! # A Genbank to GFF parser
//!
//!
//!
//! You are able to parse genbank and save as a GFF (gff3) format as well as extracting DNA sequences, gene DNA sequences (ffn) and protein fasta sequences (faa)
//!
//! You can also create new records and save as a genbank (gbk) format
//!
//! ## Detailed Explanation
//!
//!
//! The Genbank parser contains:
//!
//! BLAST records are provided by using any of these well-known programs to determine similarity
//! BLAST+ from NCBI as BLASTP, BLASTX, BLASTN, TBLASTN or TBLASTX, Diamond BLAST (protein only) or MMSeqs2
//! Several different output formats can be specified, here we provide parsers for two common output formats:
//! The BLAST+ outfmt 5 (an XML verbose format) and the outfmt 6 (the single line tabular format)
//!
//!
//!  Example to save a provided BLAST outfmt 6 single line tabular format into json
//!
//!```rust
//!    use std::io::Write;
//!    use tokio::io::BufReader as TokioBufReader;
//!    use std::io::Cursor;
//!    use microBioRust::blast::stream_outfmt6_to_json;
//!   
//!    async fn test_stream_tab_to_json() {
//!        let data = "q1  h1      99.0    10      0       0       1       10      1       10      1e-5    50";
//!        let cursor = Cursor::new(data.as_bytes());
//!        let reader = TokioBufReader::new(cursor);
//!        let res = stream_outfmt6_to_json(reader).await;
//!        println!("results are {:?}", &res);
//!    }
//!```
//! Example to create a completely new blast XML record
//!
//!```rust
//! use std::io::Write;
//! use tokio::io::BufReader as TokioBufReader;
//! use std::io::Cursor;
//! use microBioRust::blast::AsyncBlastXmlIter;
//!
//! async fn test_async_xml_iter_simple() {
//!       let xml = r#"<?xml version="1.0"?>
//!<BlastOutput>
//!  <BlastOutput_iterations>
//!    <Iteration>
//!      <Iteration_query-ID>Query_1</Iteration_query-ID>
//!      <Iteration_query-def>My query</Iteration_query-def>
//!      <Iteration_query-len>100</Iteration_query-len>
//!      <Hit>
//!        <Hit_id>gi|1</Hit_id>
//!        <Hit_def>Some hit</Hit_def>
//!        <Hit_accession>ABC123</Hit_accession>
//!        <Hit_len>100</Hit_len>
//!        <Hsp>
//!          <Hsp_bit-score>50.0</Hsp_bit-score>
//!          <Hsp_evalue>1e-5</Hsp_evalue>
//!          <Hsp_query-from>1</Hsp_query-from>
//!          <Hsp_query-to>100</Hsp_query-to>
//!          <Hsp_hit-from>1</Hsp_hit-from>
//!          <Hsp_hit-to>100</Hsp_hit-to>
//!          <Hsp_identity>90</Hsp_identity>
//!          <Hsp_align-len>100</Hsp_align-len>
//!        </Hsp>
//!      </Hit>
//!    </Iteration>
//!  </BlastOutput_iterations>
//!</BlastOutput>"#;
//!
//!        let cursor = Cursor::new(xml.as_bytes());
//!        let reader = TokioBufReader::new(cursor);
//!        let mut iter = AsyncBlastXmlIter::from_reader(reader);
//!        let next = iter.next_iteration().await;
//!
//!        assert!(next.is_some());
//!        let it = next.unwrap().unwrap();
//!        assert_eq!(it.query_id.unwrap(), "Query_1");
//!        assert_eq!(it.query_def.unwrap(), "My query");
//!        assert_eq!(it.hits.len(), 1);
//!    }
//!
//!```
#![allow(unused_imports)]
use anyhow::{Context, Result};
use async_compression::tokio::bufread::GzipDecoder as AsyncGzDecoder;
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use quick_xml::escape::unescape;
use serde::Serialize;
use serde_json::ser::Serializer as JsonSerializer;
use std::io::Cursor;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncWriteExt, BufReader};


#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct BlastTabRecord {
    pub qseqid: String,
    pub sseqid: String,
    pub pident: f32,
    pub length: u32,
    pub mismatch: Option<u32>,
    pub gapopen: Option<u32>,
    pub qstart: Option<u32>,
    pub qend: Option<u32>,
    pub sstart: Option<u32>,
    pub send: Option<u32>,
    pub evalue: f64,
    pub bitscore: f64,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Hsp {
    pub bit_score: f64,
    pub score: Option<f64>,
    pub evalue: f64,
    pub query_from: Option<u32>,
    pub query_to: Option<u32>,
    pub hit_from: Option<u32>,
    pub hit_to: Option<u32>,
    pub identity: Option<u32>,
    pub positive: Option<u32>,
    pub gaps: Option<u32>,
    pub align_len: Option<u32>,
    pub qseq: Option<String>,
    pub hseq: Option<String>,
    pub midline: Option<String>,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Hit {
    pub id: Option<String>,
    pub def: Option<String>,
    pub accession: Option<String>,
    pub len: Option<u32>,
    pub hsps: Vec<Hsp>,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct BlastXmlIteration {
    pub query_id: Option<String>,
    pub query_def: Option<String>,
    pub query_len: Option<u32>,
    pub hits: Vec<Hit>,
    pub stats: Option<Statistics>,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Statistics {
    pub db_num: Option<u32>,
    pub db_len: Option<u32>,
    pub hsp_len: Option<u32>,
    pub eff_space: Option<u32>,
    pub kappa: Option<f64>,
    pub lambda: Option<f64>,
    pub entropy: Option<u32>,
}

// first macro: For Option fields (Strings, Option<u32>, etc.)
macro_rules! read_parse_opt {
    ($self:expr, $tag:expr, $parent_opt:expr, $field:ident, $ok_ty:ty) => {{
        //read content first. This borrows 'self', but the borrow ends
        //as soon as this statement finishes.
        let res = $self.read_tag_content($tag).await;

        //now 'self' is free again and we can handle the result.
        match res {
            Ok(text) => {
                // 3. Re-borrow only the specific field we want to update.
                if let Some(parent) = &mut $parent_opt {
                    parent.$field = text.parse().ok();
                }
            }
            Err(e) => {
                return Some(std::result::Result::<$ok_ty, anyhow::Error>::Err(anyhow::Error::from(e)));
            }
        }
    }};
}
// second macro: For fields that are numerical types (e.g. bit_score, evalue)
// Uses .unwrap_or_default() to return 0.0 or 0 if parsing fails.
macro_rules! read_parse_val {
    ($self:expr, $tag:expr, $parent_opt:expr, $field:ident, $ok_ty:ty) => {{
        // 1. Read first
        let res = $self.read_tag_content($tag).await;
        
        // 2. Update second
        match res {
            Ok(text) => {
                if let Some(parent) = &mut $parent_opt {
                    parent.$field = text.parse().unwrap_or_default();
                }
            }
            Err(e) => {
                return Some(std::result::Result::<$ok_ty, anyhow::Error>::Err(anyhow::Error::from(e)));
            }
        }
    }};
}

// Async gzip-aware reader, decides on .gz suffix
pub async fn open_async_reader(path: &str) -> Result<Box<dyn AsyncBufRead + Unpin + Send>> {
    if path == "-" {
        // Stdin
        let stdin = io::stdin();
        let reader = BufReader::new(stdin);
        return Ok(Box::new(reader));
    }
    let file = tokio::fs::File::open(path).await.context("opening file")?;
    let buf = BufReader::new(file);
    if path.ends_with(".gz") {
        let dec = AsyncGzDecoder::new(buf);
        Ok(Box::new(BufReader::new(dec)))
    } else {
        Ok(Box::new(buf))
    }
}


// BLAST (-outfmt 6) format 6 single line tabular output async parser — streaming line-by-line
pub async fn stream_outfmt6_to_json<R>(reader: R) -> Result<()>
where
    R: AsyncBufRead + Unpin + Send,
{
    let mut lines = reader.lines();
    let stdout = io::stdout();
    let mut out = stdout;
    //for newline-delimited JSON writer
    while let Some(line) = lines.next_line().await? {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 3 {
            // skip incorrect line
            continue;
        }
        let rec = BlastTabRecord {
            qseqid: cols.get(0).ok_or_else(|| anyhow::anyhow!("Missing qseqid column (0)"))?.to_string(),
            sseqid: cols.get(1).ok_or_else(|| anyhow::anyhow!("Missing sseqid column (1)"))?.to_string(),
            pident: cols.get(2).ok_or_else(|| anyhow::anyhow!("Missing pident column (2)"))?.parse().unwrap_or(0.0),
            length: cols.get(3).and_then(|s| s.parse().ok()).unwrap_or(0),
            mismatch: cols.get(4).and_then(|s| s.parse().ok()),
            gapopen: cols.get(5).and_then(|s| s.parse().ok()),
            qstart: cols.get(6).and_then(|s| s.parse().ok()),
            qend: cols.get(7).and_then(|s| s.parse().ok()),
            sstart: cols.get(8).and_then(|s| s.parse().ok()),
            send: cols.get(9).and_then(|s| s.parse().ok()),
            evalue: cols.get(10).and_then(|s| s.parse().ok()).unwrap_or(0.0),
            bitscore: cols.get(11).and_then(|s| s.parse().ok()).unwrap_or(0.0),
        };
        let mut buf = Vec::new();
        serde_json::to_writer(&mut buf, &rec)?;
        buf.push(b'\n');
        tokio::io::stdout().write_all(&buf).await?;
    }
    Ok(())
}

// Async true streaming XML parser
pub struct AsyncBlastXmlIter<R>
where
    R: AsyncBufRead + Unpin,
{
    reader: Reader<R>,
    buf: Vec<u8>,
    in_iteration: bool,
    in_hit: bool,
    in_hsp: bool,
    in_stats: bool,
    cur_iteration: Option<BlastXmlIteration>,
    cur_hit: Option<Hit>,
    cur_hsp: Option<Hsp>,
    cur_stats: Option<Statistics>,
}
impl<R> AsyncBlastXmlIter<R>
where
    R: AsyncBufRead + Unpin,
{
    pub fn from_reader(r: R) -> Self {
        //let mut reader = AsyncReader::from_reader(r);
	let mut reader = Reader::from_reader(r);
        reader.config_mut().trim_text(true);
	reader.config_mut().check_end_names = false;
        Self {
            reader,
            buf: Vec::new(),
            in_iteration: false,
            in_hit: false,
            in_hsp: false,
	    in_stats: false,
            cur_iteration: None,
            cur_hit: None,
            cur_hsp: None,
	    cur_stats: None,
        }
    }
    /// Return the next iteration asynchronously (and none on EOF)
    pub async fn read_tag_content(&mut self, end_tag_name: &[u8]) -> Result<String> {
    let mut text = String::new();
    //separate buffer for the inner loop
    loop {
       match self.reader.read_event_into_async(&mut self.buf).await {
           Ok(Event::Text(e)) => {
	       //unescape XML characters
	       let raw = std::str::from_utf8(&e).context("utf8 conversion")?;
	       let escaped = unescape(raw).context("unescaping xml")?;
	       text.push_str(&escaped);
	       }
	   Ok(Event::CData(e)) => {
	       let raw = std::str::from_utf8(&e).context("utf8 conversion")?;
	       text.push_str(raw);
	       }
	   Ok(Event::End(e)) => {
	       if e.name().as_ref() == end_tag_name {
	          self.buf.clear();
		  return Ok(text);
		  }
	       }
	   Ok(Event::Eof) => {
	       return Err(anyhow::anyhow!("Unexpected EOF while reading tag content"));
	       }
	   Err(e) => return Err(anyhow::Error::from(e)),
	   _ => { () }
	   }
	   self.buf.clear();
	   }
       }
    pub async fn next_iteration(&mut self) -> Option<Result<BlastXmlIteration>> {
        loop {
            match self.reader.read_event_into_async(&mut self.buf).await {
                Ok(Event::Start(e)) => {
                    match e.name().as_ref() {
                        b"Iteration" => {
                            self.in_iteration = true;
                            self.cur_iteration = Some(BlastXmlIteration { query_id: None, query_def: None, query_len: None, hits: Vec::new(), stats: None, });
                        }
                        b"Iteration_query-def" => {
                            if self.in_iteration {
                               read_parse_opt!(self, b"Iteration_query-def", self.cur_iteration, query_def, BlastXmlIteration);
                                }
                            }
                        b"Iteration_query-ID" => {
			    if self.in_iteration {
                               read_parse_opt!(self, b"Iteration_query-ID", self.cur_iteration, query_id, BlastXmlIteration);
                                }
                            }
                        b"Iteration_query-len" => {
                            if self.in_iteration {
                               read_parse_opt!(self, b"Iteration_query-len", self.cur_iteration, query_len, BlastXmlIteration);
                                }
                            }
                        b"Hit" => {
			    //println!("DEBUG: Found Start Hit (in_iter: {})", self.in_iteration);
                            if self.in_iteration {
                                self.in_hit = true;
                                self.cur_hit = Some(Hit { id: None, def: None, accession: None, len: None, hsps: Vec::new() });
				//println!("so have created cur_hit {:?}", &self.cur_hit);
                            }
			   }
                        b"Hit_id" => {
                            if self.in_hit {
			            read_parse_opt!(self, b"Hit_id", self.cur_hit, id, BlastXmlIteration);
                                }
                            }
                        b"Hit_def" => {
                            if self.in_hit {
			            read_parse_opt!(self, b"Hit_def", self.cur_hit, def, BlastXmlIteration);
				}
			    }
                        b"Hit_accession" => {
                            if self.in_hit {
                                    read_parse_opt!(self, b"Hit_accession", self.cur_hit, accession, BlastXmlIteration);
                                }
                            }
                        b"Hit_len" => {
                            if self.in_hit {
                                     read_parse_opt!(self, b"Hit_len", self.cur_hit, len, BlastXmlIteration);
                                }
                            }
                        b"Hsp" => {
                            if self.in_hit {
                                self.in_hsp = true;
                                self.cur_hsp = Some(Hsp { bit_score: 0.0, score: None, evalue: 0.0, query_from: None, query_to: None, hit_from: None, hit_to: None, identity: None, positive: None, gaps: None, align_len: None, qseq: None, hseq: None, midline: None });
                            }
                        }
                        b"Hsp_bit-score" => {
                            if self.in_hsp {
                               read_parse_val!(self, b"Hsp_bit-score", self.cur_hsp, bit_score, BlastXmlIteration);
                               }
                            }
                        b"Hsp_score" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_score", self.cur_hsp, score, BlastXmlIteration);
                                }
                            }
                        b"Hsp_evalue" => {
                            if self.in_hsp {
                               read_parse_val!(self, b"Hsp_evalue", self.cur_hsp, evalue, BlastXmlIteration);
                                }
                            }
                        b"Hsp_query-from" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_query-from", self.cur_hsp, query_from, BlastXmlIteration);
                                }
                            }
                        b"Hsp_query-to" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_query-to", self.cur_hsp, query_to, BlastXmlIteration);
                                }
                            }
                        b"Hsp_hit-from" => {
                            if self.in_hsp {
                                read_parse_opt!(self, b"Hsp_hit-from", self.cur_hsp, hit_from, BlastXmlIteration);
                                }
                            }
                        b"Hsp_hit-to" => {
                            if self.in_hsp {
                                read_parse_opt!(self, b"Hsp_hit-to", self.cur_hsp, hit_to, BlastXmlIteration);
                                }
                            }
                        b"Hsp_identity" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_identity", self.cur_hsp, identity, BlastXmlIteration);
                                }
                            }
			b"Hsp_positive" => {
                            if self.in_hsp {              
                               read_parse_opt!(self, b"Hsp_positive", self.cur_hsp, positive, BlastXmlIteration);
                               }
			    }
                        b"Hsp_gaps" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_gaps", self.cur_hsp, gaps, BlastXmlIteration);
			       }
			    }
                        b"Hsp_qseq" => {
                            if self.in_hsp {                        
                               read_parse_opt!(self, b"Hsp_qseq", self.cur_hsp, qseq, BlastXmlIteration);
			       }
			    }
                        b"Hsp_hseq" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_hseq", self.cur_hsp, hseq, BlastXmlIteration);
			       }
			    }
                        b"Hsp_midline" => {
                            if self.in_hsp {
                               read_parse_opt!(self, b"Hsp_midline", self.cur_hsp, midline, BlastXmlIteration);
			       }
			    }
                        b"Hsp_align-len" => {
                            if self.in_hsp {
                                read_parse_opt!(self, b"Hsp_align-len", self.cur_hsp, align_len, BlastXmlIteration);
                                }
                            }
			b"Statistics" => {
			    if self.in_iteration {
			       self.in_stats = true;
			       self.cur_stats = Some(Statistics { db_num: None, db_len: None, hsp_len: None, eff_space: None, kappa: None, lambda: None, entropy: None });
			       }
			   }
                        b"Statistics_db-num" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_db-num", self.cur_stats, db_num, BlastXmlIteration);
			      }
			  }
                        b"Statistics_db-len" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_db-len", self.cur_stats, db_len, BlastXmlIteration);
			      }
			    }
                        b"Statistics_hsp-len" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_hsp-len", self.cur_stats, hsp_len, BlastXmlIteration);
			      }
			    }			    
                       b"Statistics_eff-space" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_eff-space", self.cur_stats, eff_space, BlastXmlIteration);
			      }
			    }
                       b"Statistics_kappa" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_kappa", self.cur_stats, kappa, BlastXmlIteration);
			      }
			    }
                       b"Statistics_lambda" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_lambda", self.cur_stats, lambda, BlastXmlIteration);
			      }
			    }
                       b"Statistics_entropy" => {
                            if self.in_stats {
                              read_parse_opt!(self, b"Statistics_entropy", self.cur_stats, entropy, BlastXmlIteration);
			      }
			   }
                        _ => {},
		      }
		   }
                Ok(Event::End(e)) => {
                    match e.name().as_ref() {
                        b"Hsp" => {
                            self.in_hsp = false;
                            if let Some(hsp) = self.cur_hsp.take() {
                                if let Some(hit) = &mut self.cur_hit { hit.hsps.push(hsp); }
				else { println!("DEBUG: error tried to save hsp but cur_hit is none"); }
                            }
                        }
                        b"Hit" => {
			    println!("DEBUG: found end hit");
                            self.in_hit = false;
                            if let Some(hit) = self.cur_hit.take() {
			        println!("DEBUG: Pushing hit to iteration. HSP count: {}", hit.hsps.len());
                                if let Some(iter) = &mut self.cur_iteration { iter.hits.push(hit); }
				else { println!("DEBUG: ERROR - Found </Hit>, but cur_hit is None (was it never started?)"); }
                            }
                        }
			b"Iteration_hits" => {
			    if let Some(hit) = self.cur_hit.take() {
			       if let Some(iter) = &mut self.cur_iteration {
			           iter.hits.push(hit);
				   }
			       }
			 }
                        b"Iteration" => {
			    println!("DEBUG found end iteration");
                            self.in_iteration = false;
                            if let Some(iter) = self.cur_iteration.take() {
                                return Some(Ok(iter));
                            }
                        }
			b"Statistics" => {
                            self.in_stats = false;
                            if let Some(stats) = self.cur_stats.take() {
                                if let Some(iter) = &mut self.cur_iteration {
                                    iter.stats = Some(stats);
                                    }
                                }
                        }
                        _ => {}
                    }
                },
                Ok(Event::Eof) => return None,
                Err(e) => return Some(Err(anyhow::anyhow!(e))),
                _ => {}
               }
            self.buf.clear();
	    }
        }
}


//here we use input in a format such as this in order to capture the required format of XML (5) or Tabular (6) and true or false for Json output
//see examples for further detail

//#[derive(Parser, Debug)]
//#[command(name = "blast-parsers", author, version, about = "async microBioRust BLAST parsers: for outfmt6 (single line tabular) and outfmt5 (xml)")]
//struct Cli {
//    ///Use .gz for gzip-compressed files.
//    #[arg(short, long, default_value = "-")]
//    input: String,
//    /// Format: '6' (tabular) or '5' (xml). If omitted we try to infer by file suffix only
//    #[arg(short, long)]
//    format: Option<String>,
//    /// Output newline-delimited JSON (one JSON object per record/iteration)
//    #[arg(long)]
//    json: bool,
//}

pub fn infer_format(path: &str, explicit: &Option<String>) -> String {
    if let Some(f) = explicit { return f.clone(); }
    if path.ends_with(".xml") || path.ends_with(".xml.gz") { "5".to_string() }
    else { "6".to_string() }
}

// Unit tests (async if relevant)
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tokio::io::BufReader as TokioBufReader;
    
    #[tokio::test]
    async fn test_stream_tab_to_json() {
        let data = "q1	h1	99.0	10	0	0	1	10	1	10	1e-5	50
";
        let cursor = Cursor::new(data.as_bytes());
        let reader = TokioBufReader::new(cursor);
        // run parser (writes to stdout) — here we just ensure it doesn't error
        let res = stream_outfmt6_to_json(reader).await;
        assert!(res.is_ok());
    }

    #[tokio::test]
    async fn test_async_xml_iter_simple() {
        let xml = r#"<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>My query</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Hit>
        <Hit_id>gi|1</Hit_id>
        <Hit_def>Some hit</Hit_def>
        <Hit_accession>ABC123</Hit_accession>
        <Hit_len>100</Hit_len>
        <Hsp>
          <Hsp_bit-score>50.0</Hsp_bit-score>
          <Hsp_evalue>1e-5</Hsp_evalue>
          <Hsp_query-from>1</Hsp_query-from>
          <Hsp_query-to>100</Hsp_query-to>
          <Hsp_hit-from>1</Hsp_hit-from>
          <Hsp_hit-to>100</Hsp_hit-to>
          <Hsp_identity>90</Hsp_identity>
          <Hsp_align-len>100</Hsp_align-len>
        </Hsp>
      </Hit>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"#;
        println!("here is the xml {:?}", &xml);
        let cursor = Cursor::new(xml.as_bytes());
        let reader = TokioBufReader::new(cursor);
        let mut iter = AsyncBlastXmlIter::from_reader(reader);
        let next = iter.next_iteration().await;
	println!("next is {:?}", &next);
        assert!(next.is_some());
        let it = next.unwrap().unwrap();
        assert_eq!(it.query_id.unwrap(), "Query_1");
        assert_eq!(it.query_def.unwrap(), "My query");
        assert_eq!(it.hits.len(), 1);
    }
}
