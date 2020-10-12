from typing import Dict, List
from os import path
from decimal import Decimal
from Classes.ProteinAccessionsDB import ProteinAccessionsDB
from Classes.ProteinGroup import ProteinGroup
from Classes.ReadTable import ReadTable
from Classes.Errors import ColumnNotFoundError, AccessionNotFoundError
from Classes.Sequence import Sequence


class ProteinGroupsDB(dict):
    necessaryColumns = ["Accession", "Unused"]
    sortedTableNums: List[str]
    proteinAccessionsDB: ProteinAccessionsDB
    seqDB: Dict[str, Sequence]

    def __init__(self,
                 proteinAccessionsDB: ProteinAccessionsDB,
                 seqDB: Dict[str, Sequence],
                 folder: str) -> None:
        self.proteinAccessionsDB = proteinAccessionsDB
        self.seqDB = seqDB
        self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str) -> None:
        filenames = ProteinAccessionsDB._GetProteinFilenames(folder)
        dictionary: Dict[str, Dict[str, List[str]]] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            dictionary[tableNum] = ReadTable(filename)
        self.LoadFromDict(dictionary)

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        self.sortedTableNums = sorted([*dictionary], key=lambda x: float(x))
        for tableNum, table in dictionary.items():
            self[tableNum] = []
            for columnName in self.necessaryColumns:
                if columnName not in table:
                    raise ColumnNotFoundError(
                        columnName, f"{tableNum}_ProteinSummary.txt")
            if len(table["Accession"]):
                curGroup = ProteinGroup(
                    Decimal(table["Unused"][0]), [table["Accession"][0]])
                i = 1
                while i < len(table["Accession"]):
                    if table["Accession"][i].startswith("RRRRR"):
                        break

                    if Decimal(table["Unused"][i]) != Decimal(0):
                        curGroup.accessions = sorted(curGroup.accessions)
                        self[tableNum].append(curGroup)
                        curGroup = ProteinGroup(
                            Decimal(table["Unused"][i]), [])
                    curGroup.accessions.append(table["Accession"][i])
                    i += 1
                self[tableNum].append(curGroup)
        self.CalculateRepresentatives()

    def CalculateRepresentatives(self):
        for tableNum, table in self.items():
            for group in table:
                group.representativeAccession = (
                    self.proteinAccessionsDB.GetRepresentative(
                        group.accessions, self.seqDB))

    def GetReplacementsPerTable(self) -> Dict[str, Dict[str, str]]:
        replacements: Dict[str, Dict[str, str]] = {}
        for tableNum, groups in self.items():
            replacements[tableNum] = {}
            group: ProteinGroup
            for group in groups:
                if group.representativeAccession is None:
                    raise AccessionNotFoundError(
                        "Representative accession for group"
                        f"{group.accessions} not found")
                for accession in group.accessions:
                    replacements[tableNum][accession] = (
                        group.representativeAccession)
        return replacements