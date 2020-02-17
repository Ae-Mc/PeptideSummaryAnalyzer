from typing import List, Dict, Union, Tuple
from os import listdir, path
from Classes.ReadTable import ReadTable


class ProteinTables:
    proteinReplacements: Dict[str, Dict[str, str]]
    proteinReplacementsGroups: Dict[str, Dict[str, Dict[str, int]]]
    sortedTableNums: List[str]

    def __init__(self, inputDir: str = None) -> None:

        if inputDir is not None:
            self.GetProteinSummaryReplacements(inputDir)
            self.sortedTableNums = sorted(
                [filename.split('_')[0] for filename in listdir(inputDir)
                    if filename.endswith("ProteinSummary.txt")],
                key=lambda x: float(x))
            self.RemoveReversedAccessionsFromProteinReplacements()
            self.GetProteinReplacementsGroupsPerTable()

    def GetProteinSummaryReplacements(
            self, inputDir: str):

        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]] = (
            self.GetProteinGroupsFromFiles(inputDir))
        accessions: Dict[str, Dict[str, Union[float, int]]] = (
            self.GetRepresentativeAccessionsFromGroups(groupsPerTables))
        self.proteinReplacements = (
            self.GetReplacementsPerTable(accessions, groupsPerTables))

    def RemoveReversedAccessionsFromProteinReplacements(self) -> None:

        for tableNum in self.proteinReplacements.keys():
            replacings = [key for key in self.proteinReplacements[tableNum]]
            for replacing in replacings:
                if replacing.startswith("RRRRR"):
                    del self.proteinReplacements[tableNum][replacing]

    def GetProteinGroupsFromFiles(
            self, inputDir: str) -> Dict[str, List[List[Tuple[str, float]]]]:

        groups: Dict[str, List[List[Tuple[str, float]]]] = {}
        for filename in listdir(inputDir):
            if filename.endswith("ProteinSummary.txt"):
                groups[filename.split('_')[0]] = self.GetProteinGroupsFromFile(
                    path.join(inputDir, filename))
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
                            "occurences": 0}
                    accessions[accessionName]["occurences"] += 1
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
            accessionsWithMaxUnused[group[0][0]]["occurences"])
        for accessionName, unused in sorted(group, key=lambda x: x[0]):
            curAccession = accessionsWithMaxUnused[accessionName]
            if(curAccession["unused"] >  # type: ignore
               representativeAccession[1]
               or (
                   curAccession["unused"] == representativeAccession[1] and
                   curAccession["occurences"] >  # type: ignore
                   representativeAccession[2])):
                representativeAccession = (accessionName,
                                           curAccession["unused"],
                                           curAccession["occurences"])
        return representativeAccession[0]

    def GetProteinGroupsFromFile(
            self, filename: str) -> List[List[Tuple[str, float]]]:
        fileTable = ReadTable(filename)
        return self.GetProteinGroupsFromTable(fileTable)

    def GetProteinGroupsFromTable(
            self,
            table: Dict[str, List[str]]) -> List[List[Tuple[str, float]]]:
        groups: List[List[Tuple[str, float]]] = []
        for i in range(0, len(table["Unused"])):
            if float(table["Unused"][i]) != 0:
                if len(groups):
                    groups[-1].sort()
                unused = float(table["Unused"][i])
                groups.append([])
            groups[-1].append((table["Accession"][i],
                               unused))
        if len(groups):
            groups[-1].sort()
        return groups

    def GetProteinReplacementsGroupsPerTable(self) -> None:
        self.proteinReplacementsGroups = {}
        for tableNum, replacements in self.proteinReplacements.items():
            for replaceable, replacing in replacements.items():
                if replaceable == replacing:
                    continue
                if replacing not in self.proteinReplacementsGroups:
                    self.proteinReplacementsGroups[replacing] = {}
                if(replaceable not in
                   self.proteinReplacementsGroups[replacing]):
                    self.proteinReplacementsGroups[replacing][replaceable] = {}
                    for tableNumber in self.sortedTableNums:
                        self.proteinReplacementsGroups[
                            replacing][replaceable][tableNumber] = 0
                self.proteinReplacementsGroups[
                    replacing][replaceable][tableNum] = 1
