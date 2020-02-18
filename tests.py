import unittest
from typing import Dict, List, Tuple, Any
from Classes.ProteinTables import ProteinTables


class EqualityError(Exception):

    def __init__(self, message):
        self.message = message


class Tests(unittest.TestCase):

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
        self.assertEqual(
            self.proteinTables.GetProteinGroupsFromTable(table), result)


def main():
    unittest.main()


if __name__ == "__main__":
    main()
