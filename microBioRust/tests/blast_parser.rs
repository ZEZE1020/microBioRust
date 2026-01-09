#![allow(unused_imports)]
use microBioRust::blast::{stream_outfmt6_to_json, AsyncBlastXmlIter};
use std::io::Cursor;
use tokio::io::BufReader;

/// Tests for the BLAST parsers
/// There are parsers available for two formats, -outfmt 5 (XML)
/// and -outfmt 6 (single line tabular)
/// the tests include a json format writer for a outfmt 6 line
/// and a reader and parser for XML format 5 

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
        let xml = r#"<?xml version=\"1.0\"?>
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
        let cursor = Cursor::new(xml.as_bytes());
        let reader = TokioBufReader::new(cursor);
        let mut iter = AsyncBlastXmlIter::from_reader(reader);
        let next = iter.next_iteration().await;
        assert!(next.is_some());
        let it = next.unwrap().unwrap();
        assert_eq!(it.query_id.unwrap(), "Query_1");
        assert_eq!(it.query_def.unwrap(), "My query");
        assert_eq!(it.hits.len(), 1);
    }
}
