import json
import shutil
import unittest
from pathlib import Path

import numpy as np

from pairinteraction.simulator_old.simulation import Simulation
from pairinteraction.validator.model import Model

directory = Path(__file__).parent
names = ["simulation_1", "simulation_2", "simulation_3"]


class SimulationTests(unittest.TestCase):
    def testSimulations(self):
        for name in names:
            self._oneTestSimulation(name)

    def loadModel(self, name):
        with open(f"{directory}/models/{name}.json", encoding="utf-8") as f:
            model = json.load(f)
        self.model = Model.model_validate(model)

    def _oneTestSimulation(self, name):
        self.loadModel(name)
        simulation = Simulation(self.model)
        results_list = simulation.run()
        energies_list = [results["energies"] for results in results_list]

        reference_path = f"{directory}/data/{name}__energies.txt"
        # use this for updating the reference
        # self.saveEnergiesList(reference_path, energies_list)

        reference_list = self.loadEnergiesList(reference_path)

        assert len(energies_list) == len(reference_list)
        for energies, reference in zip(energies_list, reference_list):
            assert len(energies) == len(reference)
            assert np.allclose(energies, reference, rtol=1e-4)

        self.removeCache()

    def saveEnergiesList(self, name, energies_list):
        with open(name, "w", encoding="utf-8") as f:
            for sublist in energies_list:
                f.write(" ".join(map(str, sublist)) + "\n")

    def loadEnergiesList(self, name):
        with open(name, encoding="utf-8") as f:
            return [list(map(float, line.split())) for line in f]

    def removeCache(self):
        if self.model.numerics.path_cache is not None:
            shutil.rmtree(self.model.numerics.path_cache, ignore_errors=True)
        self.model = None


if __name__ == "__main__":
    unittest.main()
