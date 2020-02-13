from os import listdir
from typing import Dict, List
from Classes.ReadTable import ReadTable
from Classes.ProteinTables import ProteinTables

""" peptideTables - словарь вида
    {
        "1.1": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение n1"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение n1"],
            ..............................................................
            "Заголовок N1": ["Значение 1", "Значение 2", ..., "Значение n1"]
        },
        "1.2": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение n2"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение n2"],
            ..............................................................
            "Заголовок N2": ["Значение 1", "Значение 2", ..., "Значение n2"]
        },
        ...,
        "I-ый файл": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение nI"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение nI"],
            ..............................................................
            "Заголовок NI": ["Значение 1", "Значение 2", ..., "Значение nI"]
        }
    }
"""


class PeptideTables:

    peptideTables: Dict[str, Dict[str, List[str]]]

    def __init__(self, inputDir=None):
        if inputDir is not None:
            self.ReadPeptideSummaries(inputDir)
            self.sortedTableNums = self.GetSortedTableNums()
            self.RemoveReversedAccessions()
            self.RemoveExcessAccessions()

    def ReadPeptideSummaries(self, inputDir: str) -> None:
        """ Считывание всех PeptideSummary файлов в словарь """
        self.peptideTables = {}
        for filename in listdir(inputDir):
            if "Peptide" in filename:
                tableNum = filename.split('_')[0]
                self.peptideTables[tableNum] = (
                    ReadTable(inputDir + '/' + filename))

    def GetSortedTableNums(self) -> List[str]:
        return sorted(self.peptideTables.keys(), key=lambda x: float(x))

    def RemoveReversedAccessions(self):

        for peptideTable in self.peptideTables.values():
            i = 0
            while i < len(peptideTable["Accessions"]):
                if peptideTable["Accessions"][i].startswith("RRRRR"):
                    break
                i += 1

            if i < len(peptideTable["Accessions"]):
                for column in peptideTable.values():
                    del column[i:]

    def RemoveExcessAccessions(self):
        for table in self.peptideTables.values():
            table["Accessions"] = [
                accession.split(';')[0] for accession in table["Accessions"]]

    def ApplyProteinReplacements(self, proteinTables: ProteinTables):

        for tableName, table in self.peptideTables.items():
            tableReplacements = proteinTables.proteinReplacements[tableName]
            for i in range(0, len(table["Accessions"])):
                if table["Accessions"][i] in tableReplacements:
                    table["Accessions"][i] = (
                        tableReplacements[table["Accessions"][i]])
