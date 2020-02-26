from typing import List, Dict, Union, Tuple
from os import listdir, path
from Classes.ReadTable import ReadTable
from Classes.Sequence import Sequence


class ProteinTables:
    proteinReplacements: Dict[str, Dict[str, str]]
    proteinReplacementsGroups: Dict[str, Dict[str, Dict[str, int]]]
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

        self.GetProteinSummaryReplacements(inputDir)
        self.GetProteinReplacementsGroupsPerTable()

    def GetProteinSummaryReplacements(
            self, inputDir: str):

        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]] = (
            self.GetProteinGroupsFromFiles(inputDir))
        representativeAccessions: Dict[str, Dict[str, Union[float, int]]] = (
            self.GetRepresentativeAccessionsFromGroups(groupsPerTables))
        self.proteinReplacements = self.GetReplacementsPerTable(
            representativeAccessions,
            groupsPerTables)

    def GetProteinGroupsFromFiles(
            self, inputDir: str) -> Dict[str, List[List[Tuple[str, float]]]]:

        groups: Dict[str, List[List[Tuple[str, float]]]] = {}
        for filename in listdir(inputDir):
            if filename.endswith("ProteinSummary.txt"):
                groups[filename.split('_')[0]] = self.GetProteinGroupsFromFile(
                    path.join(inputDir, filename))
        return groups

    def GetProteinGroupsFromFile(
            self, filename: str) -> List[List[Tuple[str, float]]]:
        fileTable = ReadTable(filename, unsafeFlag=self.unsafeReadTableFlag)
        return self.GetProteinGroupsAndRemoveReversedAccessionsFromTable(
            fileTable)

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

    def GetRepresentativeAccessionsFromGroups(
            self,
            groupsPerTables: Dict[str, List[List[Tuple[str, float]]]]
    ) -> Dict[str, Dict[str, Union[float, int]]]:
        accessions: Dict[str, Dict[str, Union[float, int]]] = {}
        for groups in groupsPerTables.values():
            for group in groups:
                for accessionName, unused in group:
                    if accessionName not in accessions:
                        accessions[accessionName] = {
                            "unused": unused,
                            "seqLen": self.seqDB[accessionName].len}
                    if unused > accessions[accessionName]["unused"]:
                        accessions[accessionName]["unused"] = unused
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
                accessionsWithMaxUnused, group)
            for accession, unused in group:
                replacements[accession] = representativeAccession
        return replacements

    def GetRepresentativeAccessionForGroup(
            self,
            accessionsWithMaxUnused: Dict[str, Dict[str, Union[float, int]]],
            group: List[Tuple[str, float]]) -> str:

        representativeAccession: Tuple[str,
                                       Union[float, int],
                                       Union[float, int]] = (
            group[0][0],
            accessionsWithMaxUnused[group[0][0]]["unused"],
            accessionsWithMaxUnused[group[0][0]]["seqLen"])
        for accessionName, unused in sorted(group, key=lambda x: x[0]):
            curAccession = accessionsWithMaxUnused[accessionName]
            if(curAccession["unused"] >  # type: ignore
               representativeAccession[1]
               or (
                   curAccession["unused"] == representativeAccession[1] and
                   curAccession["seqLen"] > representativeAccession[2])):
                representativeAccession = (accessionName,
                                           curAccession["unused"],
                                           curAccession["seqLen"])
        return representativeAccession[0]

    def GetProteinReplacementsGroupsPerTable(self) -> None:
        self.proteinReplacementsGroups = {}
        for tableNum, replacements in self.proteinReplacements.items():
            for replaceable, replacing in replacements.items():
                if replacing not in self.proteinReplacementsGroups:
                    self.proteinReplacementsGroups[replacing] = {}
                if(replaceable not in
                   self.proteinReplacementsGroups[replacing]):
                    self.proteinReplacementsGroups[replacing][replaceable] = {}
                    for tableNum in self.sortedTableNums:
                        self.proteinReplacementsGroups[
                            replacing][replaceable][tableNum] = 0
                self.proteinReplacementsGroups[
                    replacing][replaceable][tableNum] = 1
