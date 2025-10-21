from amsim import results
from pathlib import Path

from amsim import _SimulationBuilder

if __name__ == '__main__':
    result_dir = Path('/Users/karihlynsson/Github/amsimcpp/build/amsim_1610')
    builder = _SimulationBuilder()

    print(builder)
