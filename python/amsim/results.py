from pathlib import Path
from configparser import ConfigParser

from scipy import stats

import pandas as pd
import numpy as np

class SimulationResults:
    def __init__(self, output_dir: str | Path, summarise=True) -> None:
        self._output_dir: Path = Path(output_dir).expanduser().resolve()
        self._results: dict[str, pd.DataFrame] = {} 

        self._replicates = self._find_replicates()
        self._files = self._find_files(self._replicates)
        self._read_files(self._files)

        if summarise:
            self._summarise()

    def _find_replicates(self) -> list[Path]:
        return sorted([rep for rep
                            in self._output_dir.iterdir()
                            if rep.is_dir()
                           and rep.name.startswith('rep_')])

    def _find_files(self, replicates: list[Path]) -> dict[int, dict[str, Path]]:
        files = {}
        for rep in replicates:
            rep_id = int(rep.name.split('_')[1])
            files[rep_id] = { file.stem: file for file in rep.glob('*.tsv') }
        return files

    def _read_files(self, files: dict[int, dict[str, Path]]) -> None:
        for rep_id, file in files.items():
            for name, path in file.items():
                prev: pd.DataFrame = self._results.get(name, pd.DataFrame())
                cur: pd.DataFrame = pd.read_table(path)
                cur['rep'] = rep_id
                self._results[name] = pd.concat([prev, cur], ignore_index=True)

    def _summarise(self) -> None:
        self._summarised: dict[str, pd.DataFrame] = {}
        for name, data in self._results.items():
            exclude = ['rep', 'it']
            cols = data.columns.difference(exclude)

            def q025(x):
                return np.quantile(x, 0.025)

            def q975(x):
                return np.quantile(x, 0.975)

            def lci(x):
                n = len(x)
                mean = np.mean(x)
                sem = np.std(x, ddof=1) / np.sqrt(n)
                return mean - stats.t.ppf(q=0.975, df=(n - 1)) * sem

            def uci(x):
                n = len(x)
                mean = np.mean(x)
                sem = np.std(x, ddof=1) / np.sqrt(n)
                return mean + stats.t.ppf(q=0.975, df=(n - 1)) * sem

            df = (
                data
                .groupby(['it'])[cols]
                .agg(['mean', 'median', 'std', 'sem', q025, q975, lci, uci])
                .reset_index()
            )

            df.columns = ['_'.join(col).strip('_') for col in df.columns]

            df = df.melt(
                id_vars='it',
                var_name='variable',
                value_name='value'
            )

            df[['name', 'stat']] = (
                df['variable'].str.rsplit('_', n=1, expand=True)
            )

            df = (
                df.drop(columns="variable")
                  .pivot(index = ['it', 'name'], columns='stat', values='value')
            )[['mean', 'median', 'std', 'sem', 'q025', 'q975', 'lci', 'uci']]

            self._results[name] = df
