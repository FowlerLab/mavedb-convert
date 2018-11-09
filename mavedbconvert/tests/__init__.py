import os
from unittest import TestCase

import pandas as pd


__all__ = [
    'test_arg_parsing',
    'test_base',
    'test_empiric',
    'test_enrich',
    'test_enrich2',
    'test_fasta',
    'test_utilities',
    'test_filters',
    'test_validators',
    'ProgramTestCase',
]


class ProgramTestCase(TestCase):
    def setUp(self):
        self.bin = []

    def mock_multi_sheet_excel_file(self, path, data):
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        for i, di in enumerate(data):
            df = pd.DataFrame(di)
            df.to_excel(writer, sheet_name='Sheet{}'.format(i))
        writer.save()
        self.bin.append(path)

    def tearDown(self):
        for path in self.bin:
            if os.path.exists(path) and os.path.isfile(path):
                os.remove(path)
            elif os.path.exists(path) and os.path.isdir(path):
                os.removedirs(path)
