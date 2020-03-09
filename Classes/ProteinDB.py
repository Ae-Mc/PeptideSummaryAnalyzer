from os import listdir, path
from typing import Dict, List

from Classes.ProteinAccessionWithMaxUnused import ProteinAccessionWithMaxUnused
from Classes.ProteinGroup import ProteinGroup
from Classes.ReadTable import ReadTable
from Classes.Sequence import Sequence
from Classes.Errors import AccessionNotFoundError


class ProteinDB:

    inputDir: str
    proteinGroupsPerTable: Dict[str, List[ProteinGroup]]
    difficultCases: Dict[str, List[ProteinGroup]]  # Filename: [ProteinGroups]
    NecessaryColumns: List[str] = ["Unused", "Accession"]
    seqDB: Dict[str, Sequence]
    unsafeReadTableFlag: bool

    def __init__(
            self,
            inputDir: str,
            seqDB: Dict[str, Sequence],
            unsafeReadTableFlag: bool) -> None:
        self.proteinGroupsPerTable: Dict[str, List[ProteinGroup]] = {}
        self.difficultCases: Dict[str, ProteinGroup] = {}
        self.unsafeReadTableFlag = unsafeReadTableFlag
        self.seqDB = seqDB
        self.inputDir = inputDir

    def ReadDBFromFolder(self, inputDir: str = None) -> None:
        if inputDir is None:
            inputDir = self.inputDir
        for filename in listdir(inputDir):
            if filename.endswith("ProteinSummary.txt"):
                self.proteinGroupsPerTable[filename.split('_')[0]] = (
                    self._ReadDBTableFromFile(path.join(inputDir, filename)))

    def _ReadDBTableFromFile(self, filename: str) -> List[ProteinGroup]:
        rawTable = ReadTable(filename, unsafeFlag=self.unsafeReadTableFlag)
        table: List[ProteinGroup] = []
        curGroup: List[str] = []
        unused: float = 0.0

        for i in range(len(rawTable["Unused"])):
            if rawTable["Accession"][i].startswith("RRRRR"):
                break
            if float(rawTable["Unused"][i]) != 0.0:
                if len(curGroup) > 0:
                    table.append(ProteinGroup(unused, curGroup))
                unused = float(rawTable["Unused"][i])
                curGroup = []
            curGroup.append(rawTable["Accession"][i])

        if len(curGroup) > 0:
            table.append(ProteinGroup(unused, curGroup))
        return table

    def GetAccessionsBunch(self) -> Dict[str, ProteinAccessionWithMaxUnused]:
        accessionsBunch: Dict[str, ProteinAccessionWithMaxUnused] = {}

        for table in self.proteinGroupsPerTable.values():
            for group in table:
                for accession in group.accessions:
                    if (accession not in accessionsBunch):
                        # Задаём 0, чтобы неоднозначность не выставилась в
                        # True сразу
                        accessionsBunch[accession] = (
                            ProteinAccessionWithMaxUnused(accession, 0.0))
                    if accessionsBunch[accession].unused < group.unused:
                        accessionsBunch[accession].unused = group.unused
                    elif (accessionsBunch[accession].unused == group.unused):
                        accessionsBunch[accession].ambiguity = True

                    accessionsBunch[accession].occurences += 1
        return accessionsBunch

    def CalculateRepresentatives(self) -> None:
        accessionsBunch = self.GetAccessionsBunch()
        for tableNum, table in self.proteinGroupsPerTable.items():
            for group in table:
                self._CalculateRepresentativeForGroup(
                        tableNum, accessionsBunch, group)

    def _CalculateRepresentativeForGroup(
            self,
            tableNum: str,
            accessionsBunch: Dict[str, ProteinAccessionWithMaxUnused],
            group: ProteinGroup):

        candidates = self._GetCandidatesForGroup(group, accessionsBunch)
        group.representativeAccession = self._GetRepresentativeFromCandidates(
                candidates)
        self._AddGroupToDifficultCases(group, tableNum, candidates)

    def _GetCandidatesForGroup(
            self,
            group: ProteinGroup,
            accessionsBunch: Dict[str, ProteinAccessionWithMaxUnused]
            ) -> List[ProteinAccessionWithMaxUnused]:
        candidates: List[ProteinAccessionWithMaxUnused] = []
        representativeAccession = ProteinAccessionWithMaxUnused(
                group.accessions[0], 0.0)
        for accession in group.accessions:
            if (accessionsBunch[accession].unused >
                    representativeAccession.unused):
                representativeAccession = accessionsBunch[accession]
                candidates = []
            if (accessionsBunch[accession].unused ==
                    representativeAccession.unused):
                candidates.append(accessionsBunch[accession])
        return candidates

    def _GetRepresentativeFromCandidates(
            self,
            candidates: List[ProteinAccessionWithMaxUnused]) -> str:
        representativeAccession = candidates[0]
        if len(candidates) > 1:
            for candidate in candidates:
                if candidate.ambiguity:
                    return (
                        self._GetRepresentativeFromCandidatesByOccurences(
                            candidates))
            else:
                return self._GetRepresentativeFromCandidatesBySeqlens(
                        candidates)
        return representativeAccession.name

    def _GetRepresentativeFromCandidatesBySeqlens(
            self, candidates: List[ProteinAccessionWithMaxUnused]) -> str:
        representativeAccession: str = candidates[0].name
        try:
            self.seqDB[representativeAccession]
        except KeyError:
            raise AccessionNotFoundError(
                    f"Accession \"{representativeAccession}\" not found "
                    "in sequence database!")
        for candidate in candidates:
            try:
                if (self.seqDB[representativeAccession].len <
                        self.seqDB[candidate.name].len):
                    representativeAccession = candidate.name
            except KeyError:
                raise AccessionNotFoundError(
                        f"Accession \"{candidate.name}\" not found "
                        "in sequence database!")
        return representativeAccession

    def _GetRepresentativeFromCandidatesByOccurences(
            self, candidates: List[ProteinAccessionWithMaxUnused]) -> str:
        representativeAccession = candidates[0]
        for candidate in candidates:
            if candidate.occurences > representativeAccession.occurences:
                representativeAccession = candidate
        return representativeAccession.name

    def _AddGroupToDifficultCases(
            self,
            group: ProteinGroup,
            tableNum: str,
            candidates: List[ProteinAccessionWithMaxUnused]) -> None:
        if len(candidates) > 1:
            for candidate in candidates:
                if candidate.ambiguity:
                    try:
                        self.difficultCases[tableNum].append(group)
                    except KeyError:
                        self.difficultCases[tableNum] = [group]

    def GetAccessionsReplacementsPerTable(self) -> Dict[str, Dict[str, str]]:
        replacements: Dict[str, Dict[str, str]] = {}
        for tableNum, table in self.proteinGroupsPerTable.items():
            replacements[tableNum] = {}
            for group in table:
                for accession in group.accessions:
                    if group.representativeAccession is None:
                        raise AccessionNotFoundError(
                                "Representative accession for group "
                                f"{group.accessions} not found")
                    replacements[tableNum][accession] = (
                            group.representativeAccession)
        return replacements

    def GetGroupsInOutputFormat(
            self,
            groupsPerTable: Dict[str, List[ProteinGroup]]
            ) -> Dict[str, Dict[str, Dict[str, int]]]:
        """Представление групп в формате
        {
            "RID": {
                "ID1": {
                    "0.1": 0,
                    "1.1": 1,
                    .........
                    "N.n": 0,
                },
                "IDX": {
                    "0.1": 1,
                    "1.1": 0,
                    .........
                    "N.n": 0,
                }
            }
        }"""
        accessions: Dict[str, Dict[str, Dict[str, int]]] = {}
        groups: List[ProteinGroup]
        for tableNum, groups in groupsPerTable.items():
            for group in groups:
                reprAccession = group.representativeAccession
                if reprAccession is None:
                    raise AccessionNotFoundError(
                        "Representative accession for group"
                        f"{group.accessions} is not found")
                if reprAccession not in accessions:
                    accessions[reprAccession] = {}
                for accession in group.accessions:
                    if (accession not in accessions[reprAccession]):
                        accessions[reprAccession][accession] = (
                                {tableNum: 0 for tableNum in
                                    self.proteinGroupsPerTable}
                        )
                    accessions[reprAccession][accession][tableNum] = 1
        return accessions

    def GetSortedTableNums(self):
        return sorted(
                self.proteinGroupsPerTable.keys(),
                key=lambda x: float(x)
                )
