#!/usr/local/bin/python
import unittest
from io import StringIO
from typing import Dict, List, Tuple, Any
from Classes.ProteinTables import ProteinTables
from Classes.ReadTable import ReadTableFromFileObj


class ProteinTablesTest(unittest.TestCase):

    proteinTables: ProteinTables
    error: str

    def setUp(self):
        self.proteinTables = ProteinTables()

    def testGetProteinGroupsFromTable(self) -> None:
        table: Dict[str, List[str]] = {
            "Unused": [
                "10",    "0",     "0",     "2",     "0",     "0"
            ],
            "Accession": [
                "Acc01", "Acc02", "Acc03", "Acc06", "Acc05", "Acc04"
            ]}

        result: List[List[Tuple[Any, ...]]] = [
            [("Acc01", 10.0), ("Acc02", 10.0), ("Acc03", 10.0)],
            [("Acc04", 02.0), ("Acc05", 02.0), ("Acc06", 02.0)]
        ]
        self.assertListEqual(
            self.proteinTables.GetProteinGroupsFromTable(table), result)


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
