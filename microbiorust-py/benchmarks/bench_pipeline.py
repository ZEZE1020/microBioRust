import os
import time
import subprocess
from Bio import SeqIO
from microbiorust import gbk
from codecarbon import OfflineEmissionsTracker

class PipelineSuite:
    """
    Benchmarks for microbiorust-py vs BioPython.
    Measures:
      - Parsing Time
      - Peak Memory
      - Latency
      - Energy Consumption (CodeCarbon)
    """

    timeout = 300  # 5 minutes per benchmark

    # Parameterizing creates 4 lines on one graph: Rust/Python × Interactive/Pipeline
    params = [['rust', 'python'], ['interactive', 'pipeline']]
    param_names = ['engine', 'context']

    def setup(self, engine, context):
        bench_dir = os.path.dirname(__file__)
        self.filepath = os.path.join(bench_dir, "Rhiz3841.gbk.gb")

        # CLI scripts for pipeline mode
        self.rust_cli = os.path.join(bench_dir, "rust_via_python_countgbk2faa.py")
        self.python_cli = os.path.join(bench_dir, "bp_equivalent_gbktofaacount.py")

        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Missing benchmark file: {self.filepath}")

        # Initialize energy dictionary if not already present
        if not hasattr(self, "_energy_joules"):
            self._energy_joules = {}

    # --- CORE LOGIC WITH CODECARBON ---
    def _run_logic(self, engine, context):
        """
        Routes execution based on engine/context and tracks energy with CodeCarbon.
        Stores last measured energy per engine in self._energy_joules.
        """
        tracker = OfflineEmissionsTracker(measure_power_secs=1, log_level="CRITICAL", country_iso_code="USA", constant_carbon_intensity=475)
        tracker.start()

        try:
            # --- DISPATCH ---
            if context == 'interactive':
                if engine == 'rust':
                    result = gbk.gbk_to_faa_count(self.filepath)
                else:
                    count = 0
                    for record in SeqIO.parse(self.filepath, "genbank"):
                        for feature in record.features:
                            if feature.type == "CDS":
                                _ = feature.extract(record.seq).translate()
                                count += 1
                    result = count

            elif context == 'pipeline':
                if engine == 'rust':
                    result = subprocess.run(
                        ["python", self.rust_cli, self.filepath],
                        capture_output=True, check=True, timeout=self.timeout
                    ).returncode
                else:
                    result = subprocess.run(
                        ["python", self.python_cli, self.filepath],
                        capture_output=True, check=True, timeout=self.timeout
                    ).returncode
            else:
                raise ValueError(f"unknown context: {context}")

        finally:
            tracker.stop()
            energy_kwh = tracker.final_emissions_data.energy_consumed
            # Store energy per engine in Joules
            self._energy_joules[f"{engine}_{context}"] = energy_kwh * 3_600_000 # Joules

        return result

    # --- 1. PRIMARY TIME BENCHMARK (ASV automatic) ---
    def time_process_all(self, engine, context):
        """Measures parsing speed (seconds)."""
        self._run_logic(engine, context)

    # --- 2. PEAK MEMORY BENCHMARK ---
    def peakmem_process_all(self, engine, context):
        """Measures maximum resident set size (RAM)."""
        self._run_logic(engine, context)

    # --- 3. PARSING LATENCY (single-pass) ---
    def track_parsing_latency(self, engine, context):
        """Records single-pass parsing latency (seconds)."""
        start = time.perf_counter()
        self._run_logic(engine, context)
        end = time.perf_counter()
        return end - start
    track_parsing_latency.unit = "s"

    # --- 4. ENERGY BENCHMARK ---
    def track_energy(self, engine, context):
        """
        Returns energy consumed (Joules) for this engine/context using CodeCarbon.
        ASV can plot Rust vs Python separately.
        """
        self._run_logic(engine, context)
        return self._energy_joules[f"{engine}_{context}"]
    track_energy.unit = "J"

