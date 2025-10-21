from os import listdir
from pathlib import Path

import numpy as np
import pandas as pd


class SimulationResults:
    def __init__(self, result_dir: Path | str) -> None:
        self.result_dir: Path = Path(result_dir).resolve()
        self.files: list[Path] = [Path(result_dir, f) for f in listdir(self.result_dir)]
        self.data: dict[Path, pd.DataFrame] = {}

    def load_results(self):
        for f in self.files:
            self.data[f] = pd.read_table(f, sep="\t")
        print(self.data.values())
        pass

    def summarise_results(self):
        pass
