# https://johnnado.com/pyqt-qtest-example/
# https://github.com/jmcgeheeiv/pyqttestexample
import json
import os
import sys
import unittest
import zipfile

import numpy as np
import scipy.io
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication

import pairinteraction_gui.pairinteraction.app as piGui

app = QApplication(sys.argv)


class PairinteractionGuiTest(unittest.TestCase):
    def setUp(self):
        self.form = piGui.MainWindow()
        self.form.ui.action_sconf_reset.trigger()
        self.form.ui.action_pconf_reset.trigger()

    def testFieldCalcButton(self):
        self.form.ui.spinbox_system_cores.setValue(1)
        for x in "xyz":
            for minmax in ["min", "max"]:
                getattr(self.form.ui, f"lineedit_system_{minmax}E{x}").setText("1")
                getattr(self.form.ui, f"lineedit_system_{minmax}B{x}").setText("1.5")
        self.form.ui.spinbox_system_steps.setValue(1)

        self._testEnergies(0, "Field", dE=3)

    def testPotentialCalcButton(self):
        self.form.ui.lineedit_system_minR.setText("20")
        self.form.ui.lineedit_system_maxR.setText("5")
        self.form.ui.spinbox_system_steps.setValue(2)

        self._testEnergies(2, "Potential", dE=0.3)

    def _testEnergies(self, idx, ref_data, dE, dE_tol=1e-3, use_python_api="both"):
        if use_python_api == "both":
            for use_python_api in [False, True]:
                self._testEnergies(idx, ref_data, dE, use_python_api=use_python_api)
            return
        self.form.ui.checkbox_use_python_api.setChecked(use_python_api)
        self.form.autosetSymmetrization()

        if idx == 0:
            widget_calc = self.form.ui.pushbutton_field1_calc
            widget_save = self.form.ui.pushbutton_field1_save
        elif idx == 2:
            widget_calc = self.form.ui.pushbutton_potential_calc
            widget_save = self.form.ui.pushbutton_potential_save

        # Run simulation
        QTest.mouseClick(widget_calc, Qt.LeftButton)
        while self.form.timer.isActive():
            self.form.checkForData()

        # Save current data
        path = "reference_data/"
        self.form.forceFilename = path + "tmp"
        QTest.mouseClick(widget_save, Qt.LeftButton)

        # Load current and reference data
        data = {}
        sconfig = {}
        for k in ["tmp", ref_data]:
            with zipfile.ZipFile(path + k, "r") as zip_file:
                with zip_file.open("data.mat") as f:
                    data[k] = scipy.io.loadmat(f)
                with zip_file.open("settings.sconf") as f:
                    sconfig[k] = json.load(f)
        os.remove(path + "tmp")

        # Check if configs match
        for k, v in sconfig[ref_data].items():
            assert sconfig["tmp"][k] == v

        # Check if central eigenvalues (+/- dE) match
        for i in range(len(data[ref_data]["eigenvalues"])):
            Es = {k: np.array(mat["eigenvalues"])[i] for k, mat in data.items()}
            Es = {k: E[np.abs(E) < dE] for k, E in Es.items()}
            diff_rel = np.abs(Es[ref_data] - Es["tmp"])
            assert np.all(diff_rel <= dE_tol)

    def tearDown(self):
        # Calculation runs in the background. Wait for it to finish.
        if self.form.thread.isRunning():
            self.form.thread.wait()
        # Close any pipes and wait for subprocess to exit.
        if self.form.proc:
            self.form.proc.stdout.close()
            self.form.proc.wait()


if __name__ == "__main__":
    unittest.main()
