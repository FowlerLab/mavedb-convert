import os
import shutil
from unittest import TestCase
from tempfile import TemporaryDirectory

import pandas as pd


__all__ = [
    "test_base",
    "test_empiric",
    "test_enrich",
    "test_enrich2",
    "test_fasta",
    "test_utilities",
    "test_filters",
    "test_validators",
    "ProgramTestCase",
]


# TODO: think up a better name for this class
# TODO: remove the old self.bin stuff
class ProgramTestCase(TestCase):
    def setUp(self):
        self._data_dir = TemporaryDirectory()  # store the object
        self.data_dir = os.path.join(
            self._data_dir.name, "data"
        )  # store the directory path
        shutil.copytree(
            src=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"),
            dst=self.data_dir,
        )
        self.bin = []

    def mock_multi_sheet_excel_file(self, path, data):
        writer = pd.ExcelWriter(path, engine="xlsxwriter")
        for i, di in enumerate(data):
            df = pd.DataFrame(di)
            df.to_excel(writer, sheet_name="Sheet{}".format(i), index=False)
        writer.save()
        self.bin.append(path)

    def tearDown(self):
        self._data_dir.cleanup()
        for path in self.bin:
            if os.path.exists(path) and os.path.isfile(path):
                os.remove(path)
            elif os.path.exists(path) and os.path.isdir(path):
                os.removedirs(path)
