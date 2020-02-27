#!/bin/env python3
import unittest
from io import StringIO
from typing import List, Tuple, Any
from Classes.ProteinTables import ProteinTables
from Classes.ReadTable import ReadTableFromFileObj


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

        result: List[List[Tuple[Any, ...]]] = {
            "0.1": [
                ["Acc01", "Acc02", "Acc03"],
                ["Acc04", "Acc05", "Acc06"]
            ]}
        self.assertDictEqual(
            self.proteinTables.GetPureGroupsListsFromRawTables(),
            result)


class ReadTableTest(unittest.TestCase):

    def testSafeReadTableException(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\n" +
            "5\t6\t7\t8\t9\n" +
            "a\tb\tc\td\n")
        with self.assertRaises(IndexError):
            ReadTableFromFileObj(inputFile, '\t', False)

    def testUnsafeReadTable(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\t\n" +
            "5\t6\t7\t8\t\n" +
            "a\tb\tc\td\n")
        self.assertDictEqual(
            ReadTableFromFileObj(inputFile, '\t', True),
            {
                'A': ["1", "5", "a"],
                'B': ["2", "6", "b"],
                'C': ["3", "7", "c"],
                'D': ["4", "8", "d"]
            })

    def testSafeReadTable(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\n" +
            "5\t6\t7\t8\n" +
            "a\tb\tc\td\n")
        self.assertDictEqual(
            ReadTableFromFileObj(inputFile, '\t', False),
            {
                'A': ["1", "5", "a"],
                'B': ["2", "6", "b"],
                'C': ["3", "7", "c"],
                'D': ["4", "8", "d"]
            })


def main():
    unittest.main()


if __name__ == "__main__":
    main()
