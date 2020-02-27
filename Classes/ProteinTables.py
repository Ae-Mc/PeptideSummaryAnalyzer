from typing import List, Dict, Union, Tuple
from os import listdir, path
from Classes.ReadTable import ReadTable
from Classes.Sequence import Sequence


class ProteinTables:
    rawProteinTables: Dict[str, Dict[str, List[str]]]
    proteinReplacements: Dict[str, Dict[str, str]]
    proteinReplacementsGroups: Dict[str, Dict[str, Dict[str, int]]]
    interimGroupsLists: Dict[str, List[List[str]]]
    sortedTableNums: List[str]
    unsafeReadTableFlag: bool
    seqDB: Dict[str, Sequence]

    def __init__(self,
                 inputDir: str,
                 seqDB: Dict[str, Sequence],
                 unsafeReadTableFlag: bool) -> None:

        self.unsafeReadTableFlag = unsafeReadTableFlag
        self.seqDB = seqDB

        self.sortedTableNums = sorted(
            [filename.split('_')[0] for filename in listdir(inputDir)
                if filename.endswith("ProteinSummary.txt")],
            key=lambda x: float(x))

        self.GetRawPeptideTables(inputDir)
        self.GetProteinSummaryReplacements()
        self.GetProteinReplacementsGroupsPerTable()

    def GetRawPeptideTables(self, inputDir: str):
        self.rawProteinTables = {}
        for filename in listdir(inputDir):
            if filename.endswith("ProteinSummary.txt"):
                tableNum = filename.split('_')[0]
                self.rawProteinTables[tableNum] = ReadTable(
                    path.join(inputDir, filename),
                    unsafeFlag=self.unsafeReadTableFlag)

    def GetProteinSummaryReplacements(
            self):
        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]] = (
            self.GetProteinGroupsFromFiles())
        accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]] = (
            self.GetAccessionsMaxUnusedBunch())
        self.proteinReplacements = self.GetReplacementsPerTable(
            accessionsWithMaxUnused,
            groupsPerTables)

    def GetProteinGroupsFromFiles(
            self) -> Dict[str, List[List[Tuple[str, float]]]]:
        groups: Dict[str, List[List[Tuple[str, float]]]] = {}
        for tableNum, table in self.rawProteinTables.items():
            groups[tableNum] = (
                self.GetProteinGroupsAndRemoveReversedAccessionsFromTable(
                    table))
        return groups

    def GetProteinGroupsAndRemoveReversedAccessionsFromTable(
            self,
            table: Dict[str, List[str]]) -> List[List[Tuple[str, float]]]:
        groups: List[List[Tuple[str, float]]] = []
        for i in range(0, len(table["Unused"])):
            if table["Accession"][i].startswith("RRRRR"):
                if len(groups):
                    groups[-1].sort()
                break

            if float(table["Unused"][i]) != 0:
                if len(groups):
                    groups[-1].sort()
                unused = float(table["Unused"][i])
                groups.append([])
            groups[-1].append((table["Accession"][i], unused))

        if len(groups):
            groups[-1].sort()
        return groups

    def GetAccessionsMaxUnusedBunch(
            self) -> Dict[str, Dict[str, Union[float, int]]]:
        accessions: Dict[str, Dict[str, Union[float, int]]] = {}
        for table in self.rawProteinTables.values():
            i = 0
            while i < len(table["Unused"]):
                if float(table["Unused"][i]) != 0.0:
                    unused = float(table["Unused"][i])

                accessionName = table["Accession"][i]
                if accessionName.startswith("RRRRR"):
                    break
                if accessionName not in accessions:
                    accessions[accessionName] = {
                        "unused": unused,
                        "seqLen": self.seqDB[accessionName].len}
                if unused > accessions[accessionName]["unused"]:
                    accessions[accessionName]["unused"] = unused
                i += 1
        return accessions

    def GetReplacementsPerTable(
            self,
            accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]],
            groupsPerTables: Dict[str, List[List[Tuple[str, float]]]]
    ) -> Dict[str, Dict[str, str]]:
        replacementsPerTable: Dict[str, Dict[str, str]] = {}
        for tableName, groups in groupsPerTables.items():
            replacementsPerTable[tableName] = self.GetReplacementsForGroups(
                accessionsWithMaxUnused, groups)
        return replacementsPerTable

    def GetReplacementsForGroups(
            self,
            accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]],
            groups: List[List[Tuple[str, float]]]
    ) -> Dict[str, str]:

        replacements = {}
        for group in groups:
            representativeAccession = self.GetRepresentativeAccessionForGroup(
                accessionsWithMaxUnused,
                list(map(lambda x: x[0], group)))
            for accession, unused in group:
                replacements[accession] = representativeAccession
        return replacements

    def GetRepresentativeAccessionForGroup(
            self,
            accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]],
            group: List[str]) -> str:
        representativeAccession: Tuple[str,
                                       Union[float, int],
                                       Union[float, int]] = (
            group[0],
            accessionsWithMaxUnused[group[0]]["unused"],
            accessionsWithMaxUnused[group[0]]["seqLen"])
        for accessionName in sorted(group):
            curAccession = accessionsWithMaxUnused[accessionName]
            if(curAccession["unused"] > representativeAccession[1] or (
                   curAccession["unused"] == representativeAccession[1] and
                   curAccession["seqLen"] > representativeAccession[2])):
                representativeAccession = (accessionName,
                                           curAccession["unused"],
                                           curAccession["seqLen"])
        return representativeAccession[0]

    def GetProteinReplacementsGroupsPerTable(self) -> None:
        groupsPerTables: Dict[str, List[List[str]]] = (
            self.GetPureGroupsListsFromRawTables())
        accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]] = (
            self.GetAccessionsMaxUnusedBunch())
        groupsPerTablesWithRepresentativeAccessions = (
            self.GetGroupsPerTablesWithRepresentativeAccessions(
                accessionsWithMaxUnused, groupsPerTables))
        accessionsPresenceBunch = self.GetAccessionsPresenceBunch()
        self.proteinReplacementsGroups = (
            self.ConvertGroupsPerTableToTablesPerAccession(
                groupsPerTablesWithRepresentativeAccessions,
                accessionsPresenceBunch))

    def GetPureGroupsListsFromRawTables(self) -> Dict[str, List[List[str]]]:
        groupsPerTables: Dict[str, List[List[str]]] = {}
        for tableNum, table in self.rawProteinTables.items():
            groupsPerTables[tableNum] = []
            for i in range(len(table["Unused"])):
                if table["Accession"][i].startswith("RRRRR"):
                    break
                if float(table["Unused"][i]) != 0.0:
                    groupsPerTables[tableNum].append([])
                    if len(groupsPerTables[tableNum]):
                        groupsPerTables[tableNum][-1].sort()

                groupsPerTables[tableNum][-1].append(table["Accession"][i])
                i += 1
            if len(groupsPerTables[tableNum]):
                groupsPerTables[tableNum][-1].sort()
        return groupsPerTables

    def GetGroupsPerTablesWithRepresentativeAccessions(
            self,
            accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]],
            groupsPerTables: Dict[str, List[List[str]]]
    ) -> Dict[str, Dict[str, List[str]]]:
        groupsPerTablesWithRepresentativeAccessions: Dict[
            str, Dict[str, List[str]]] = {}
        for tableNum, groups in groupsPerTables.items():
            groupsPerTablesWithRepresentativeAccessions[tableNum] = {}
            curTable = groupsPerTablesWithRepresentativeAccessions[tableNum]
            for group in groups:
                curTable[self.GetRepresentativeAccessionForGroup(
                    accessionsWithMaxUnused, group)] = group
        return groupsPerTablesWithRepresentativeAccessions

    def GetAccessionsPresenceBunch(self) -> Dict[str, Dict[str, int]]:
        accessionsPresenceBunch: Dict[str, Dict[str, int]] = {}
        for tableNum, table in self.rawProteinTables.items():
            for i in range(len(table["Accession"])):
                accessionName = table["Accession"][i]
                if accessionName not in accessionsPresenceBunch:
                    accessionsPresenceBunch[accessionName] = {
                        tableNumber: 0 for tableNumber in self.sortedTableNums
                    }
                accessionsPresenceBunch[accessionName][tableNum] = 1
        return accessionsPresenceBunch

    def ConvertGroupsPerTableToTablesPerAccession(
            self,
            groupsPerTablesWithRepresentativeAccessions: Dict[
                str, Dict[str, List[str]]],
            accessionsPresenceBunch: Dict[str, Dict[str, int]]
    ) -> Dict[str, Dict[str, Dict[str, int]]]:
        groupsWithRepresentativeAccession: Dict[
            str, Dict[str, Dict[str, int]]] = {}
        for tableNum, groups in (
                groupsPerTablesWithRepresentativeAccessions.items()):
            for representativeAccession, group in groups.items():
                groupsWithRepresentativeAccession[representativeAccession] = {}
                curAccessionGroup = (
                    groupsWithRepresentativeAccession[representativeAccession])
                for accession in group:
                    curAccessionGroup[accession] = (
                        accessionsPresenceBunch[accession])
        return groupsWithRepresentativeAccession
