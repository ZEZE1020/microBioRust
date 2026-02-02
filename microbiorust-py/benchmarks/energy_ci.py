#!/usr/bin/python

import json
from pathlib import Path

class EnergyRustCI:
    unit = "J"
    version = 1

    def track_energy(self):
        path = Path("eco-ci-rust.json")
        if not path.exists():
            return None  # ASV will mark as missing locally
        return json.loads(path.read_text())["energy_joules"]


class EnergyPythonCI:
    unit = "J"
    version = 1

    def track_energy(self):
        path = Path("eco-ci-python.json")
        if not path.exists():
            return None
        return json.loads(path.read_text())["energy_joules"]

