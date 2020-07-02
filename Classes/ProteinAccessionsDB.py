from typing import Dict, List
from os import listdir, path
from Classes.ReadTable import ReadTable
from Classes.Errors import ColumnNotFoundError
from Classes.ProteinAccession import ProteinAccession
from decimal import Decimal


class ProteinAccessionsDB(dict):
    necessaryColumns = ["Accession", "N", "Unused"]

    def __init__(self, folder: str = None):
        if folder is not None:
            self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str):
        filenames = self._GetProteinFilenames(folder)
        dictionary: Dict[str, Dict[str, List[str]]] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            dictionary[tableNum] = ReadTable(filename)
        self.LoadFromDict(dictionary)

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        for tableNum, table in dictionary.items():
            for columnName in self.necessaryColumns:
                if columnName not in table:
                    raise ColumnNotFoundError(columnName,
                                              f"{tableNum}_ProteinSummary.txt")
            curUnused: Decimal = Decimal(0)
            i = 0
            while i < len(table["Accession"]):
                if table["Accession"][i].startswith("RRRRR"):
                    break

                if Decimal(table["Unused"][i]) != Decimal(0):
                    curUnused = Decimal(table["Unused"][i])
                accession = ProteinAccession(table["Accession"][i],
                                             curUnused)
                if accession.name in self:
                    self[accession.name].unused = max(
                        self[accession.name].unused, accession.unused)
                    self[accession.name].occurences += 1
                else:
                    self[accession.name] = accession
                i += 1

    def GetRepresentative(self, accessions: List[str]) -> str:
        representativeAccession = self[accessions[0]]
        for accession in accessions:
            if(
                self[accession].unused > representativeAccession.unused or
                (self[accession].unused == representativeAccession.unused and
                 self[accession].occurences >
                 representativeAccession.occurences)):
                representativeAccession = self[accession]
        return representativeAccession.name

    @staticmethod
    def _GetProteinFilenames(folder: str) -> List[str]:
        filenames: List[str] = []
        for filename in listdir(folder):
            if filename.endswith("ProteinSummary.txt"):
                filenames.append(path.join(folder, filename))
        return filenames
