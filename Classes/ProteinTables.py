from typing import List, Dict, Union, Tuple
from Classes.ReadTable import ReadTable


class ProteinTables:
    proteinReplacements: Dict[str, Dict[str, str]]
    proteinReplacementsGroups: Dict[str, Dict[str, Dict[str, int]]]

    def __init__(self, inputDir=None):
        if inputDir is not None:
            self.GetProteinSummaryReplacements(inputDir)
            self.RemoveReversedAccessionsFromProteinReplacements()

    def GetProteinSummaryReplacements(
            self, inputDir: str):

        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]] = (
            self.GetProteinGroupsFromFiles(inputDir))
        accessions: Dict[str, Dict[str, Union[float, int]]] = (
            self.GetAccessionsWithMaxUnusedFromProteinGroups(groupsPerTables))
        self.proteinReplacements = (
            GetReplacementsPerTable(accessions, groupsPerTables))

    def GetProteinGroupsFromFiles(
            self, inputDir: str) -> Dict[str, List[List[Tuple[str, float]]]]:

        groups: Dict[str, List[List[Tuple[str, float]]]] = {}
        for filename in listdir(inputDir):
            if filename.endswith("ProteinSummary.txt"):
                groups[filename.split('_')[0]] = self.GetProteinGroupsFromFile(
                    inputDir + '/' + filename)
        return groups

    def GetAccessionsWithMaxUnusedFromProteinGroups(
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
            replacementsPerTable[tableName] = GetReplacementsForGroups(
                accessionsWithMaxUnused, groups)
        return replacementsPerTable

    def GetProteinGroupsFromFile(
            self, filename: str) -> List[List[Tuple[str, float]]]:
        fileTable = ReadTable(filename)
        groups: List[List[Tuple[str, float]]] = []
        for i in range(0, len(fileTable["N"])):
            if float(fileTable["Unused"][i]) != 0:
                if len(groups):
                    groups[-1].sort()
                unused = float(fileTable["Unused"][i])
                groups.append([])
            groups[-1].append((fileTable["Accession"][i],
                               unused))
        return groups

