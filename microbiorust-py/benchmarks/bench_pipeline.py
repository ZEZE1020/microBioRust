import os
import time
import subprocess
from Bio import SeqIO
import microbiorust

class PipelineSuite:
    """
    Benchmarks for microbiorust-py vs BioPython.
    Measures performance across Time, Peak RAM, and Latency.
    """

    timeout = 300  # 5 minute timeout per benchmark

    # Parameterizing creates 4 lines on one graph for comparison
    params = [['rust', 'python'], ['interactive','pipeline']]
    param_names = ['engine', 'context']

    def setup(self, engine, context):
        bench_dir = os.path.dirname(__file__)
        self.filepath = os.path.join(bench_dir, "Rhiz3841.gbk.gb")
        
        # CLI Paths
        self.rust_cli = os.path.join(bench_dir, "rust_via_python_countgbk2faa.py") # Assuming installed/in path
        self.python_cli = os.path.join(bench_dir, "bp_equivalent_gbktofaacount.py")

        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Missing benchmark file: {self.filepath}")

    # --- 1. PRIMARY TIME BENCHMARK (Automatic) ---
    def time_process_all(self, engine, context):
        """Measures parsing speed (Standard ASV timing)."""
        self._run_logic(engine, context)

    # --- 2. PEAK MEMORY BENCHMARK ---
    def peakmem_process_all(self, engine, context):
        """Measures maximum resident set size (RAM)."""
        self._run_logic(engine, context)


    # --- 4. PARSING LATENCY (Single Pass) ---
    def track_parsing_latency(self, engine, context):
        """
        Records single-pass parsing latency. 
        Useful for trending analysis over git commits.
        """
        start = time.perf_counter()
        self._run_logic(engine, context)
        end = time.perf_counter()
        return end - start
    
    track_parsing_latency.unit = "s"

    # --- CORE LOGIC DISPATCHER ---
    def _run_logic(self, engine, context):
        """Routes execution based on the current benchmark parameter."""
        if context == 'interactive':
            if engine == 'rust':
                return microbiorust.gbk_to_faa_count(self.filepath)
            else:
                # Industry standard BioPython streaming approach
                count = 0
                for record in SeqIO.parse(self.filepath, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            # Perform translation to ensure fair CPU load
                            _ = feature.extract(record.seq).translate()
                            count += 1
                return count

        elif context == 'pipeline':
            if engine == 'rust':
                result = subprocess.run(["python", self.rust_cli, self.filepath], capture_output=True, check=True, timeout=self.timeout)
                return result.returncode
            else:
                result = subprocess.run(["python", self.python_cli, self.filepath], capture_output=True, check=True, timeout=self.timeout)
                return result.returncode
        else:
            raise ValueError(f"unknown context: {context}")
