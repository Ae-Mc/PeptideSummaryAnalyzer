from typing import Dict, List, Tuple, Any
import json
from Classes.ProteinTables import ProteinTables


class EqualityError(Exception):

    def __init__(self, message):
        self.message = message


class Tests:

    proteinTables: ProteinTables
    error: str

    def __init__(
            self,
            proteinTables: ProteinTables) -> None:
        self.proteinTables = proteinTables

    def TestGetProteinGroupsFromTable(
            self,
            table: Dict[str, List[str]] = None,
            result: List[List[Tuple[str, float]]] = None) -> None:
        if table is None:
            table = {
                "Unused": [
                    "10",    "0",     "0",     "2",     "0",     "0"
                ],
                "Accession": [
                    "Acc01", "Acc02", "Acc03", "Acc06", "Acc05", "Acc04"
                ]}
        if result is None:
            result = [
                [("Acc01", 10.0), ("Acc02", 10.0), ("Acc03", 10.0)],
                [("Acc04", 02.0, 2), ("Acc05", 02.0), ("Acc06", 02.0)]
            ]
        self.CompareLists(
            self.proteinTables.GetProteinGroupsFromTable(table),
            result)

    @staticmethod
    def CompareDicts(dict1: Dict[Any, Any],
                     dict2: Dict[Any, Any]) -> bool:
        return False

    @staticmethod
    def CompareLists(list1: List[Any],
                     list2: List[Any]) -> bool:
        if len(list1) != len(list2):
            raise(EqualityError(
                f"Lists lengthes not equal:\n\t{list1}\n\t{list2}"))
        i = 0
        for i in range(0, len(list1)):
            if not isinstance(list1[i], type(list2[i])):
                raise(EqualityError(
                    f"Lists elements types on index {i} not equal:\n " +
                    f"{type(list1[i])}\n {type(list2[i])}"))
            if isinstance(list1[i], (tuple, list)):
                Tests.CompareLists(list1[i], list2[i])
            elif isinstance(list1[i], dict):
                Tests.CompareDicts(list1[i], list2[i])
            elif list1[i] != list2[i]:
                raise(EqualityError(
                    f"Lists elements on index {i} not equal:\n " +
                    f"{list1[i]} != {list2[i]}"))


def main():
    tests = Tests(ProteinTables("Input/"))
    try:
        tests.TestGetProteinGroupsFromTable()
    except EqualityError as e:
        print(e.message)


if __name__ == "__main__":
    main()
