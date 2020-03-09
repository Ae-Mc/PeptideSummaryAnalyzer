import unittest
from typing import Dict, List
from Classes.ProteinTables import ProteinTables


class ProteinTablesTest(unittest.TestCase):

    proteinTables: ProteinTables
    error: str

    def setUp(self):
        self.proteinTables = ProteinTables(None, None, None)

    def testGetProteinGroupsFromTable(self) -> None:
        self.proteinTables.rawProteinTables["0.1"] = {
            "Unused": [
                "10",    "0",     "0",     "2",     "0",     "0"
            ],
            "Accession": [
                "Acc01", "Acc02", "Acc03", "Acc06", "Acc05", "Acc04"
            ]}

        result: Dict[str, List[List[str]]] = {
            "0.1": [
                ["Acc01", "Acc02", "Acc03"],
                ["Acc04", "Acc05", "Acc06"]
            ]}
        self.assertDictEqual(
            self.proteinTables.GetPureGroupsListsFromRawTables(),
            result)
