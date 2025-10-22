from os import listdir
from pathlib import Path

import numpy as np
import pandas as pd


class SimulationResults:
    def __init__(self, output_dir: str | Path) -> None:
        self._output_dir = Path(output_dir)
        print(output_dir)
