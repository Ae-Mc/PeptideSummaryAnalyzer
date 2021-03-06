from decimal import localcontext, Decimal
from typing import Dict, Tuple, List, Optional, Set
from os import makedirs, path
from Classes.AccessionTables import AccessionTables
from Classes.Sequence import Sequence
from Classes.Accession import Accession
from Classes.ProteinGroupsDB import ProteinGroupsDB
from Classes.ProteinGroup import ProteinGroup
from Classes.Errors import RepresentativeAccessionNotFoundError
from Classes.Input import Input


class Output:
    """Класс для вывода результатов в файлы

    Attributes:
        outputDirPath: путь для выходных файлов
        seqDB: база данных последовательностей Accession
        accessionTables: таблицы с Accession вида: {
                "Номер таблицы": {
                    "Имя Accession": Accession
                }
            }
        proteinGroupsDB: база данных с Protein группами
    """
    inputParams: Input
    _accessionsBunch: Dict[str, Dict[str, Accession]]
    seqDB: Dict[str, Sequence]
    accessionTables: AccessionTables
    proteinGroupsDB: Optional[ProteinGroupsDB]
    formattedProteinGroups: Optional[
        Dict[str, Dict[str, List[Decimal]]]] = None

    def __init__(
            self,
            inputParams: Input,
            seqDB: Dict[str, Sequence],
            accessionTables: AccessionTables,
            proteinGroupsDB: ProteinGroupsDB = None) -> None:

        self.seqDB: Dict[str, Sequence] = seqDB
        self.inputParams: str = inputParams
        self.accessionTables: AccessionTables = accessionTables
        self.proteinGroupsDB = proteinGroupsDB
        if proteinGroupsDB:
            self.ConvertProteinGroupsToOutputFormat()
        self.GenerateOutputFiles()

    def ConvertProteinGroupsToOutputFormat(self) -> None:
        """Представление групп в виде, удобном для вывода {
                "RID": {
                    "ID1": {
                        "0.1": 0,
                        "1.1": RID_1.1_Unused,
                        .........
                        "N.n": 0,
                    },
                    "IDX": {
                        "0.1": RID_0.1_Unused,
                        "1.1": 0,
                        .........
                        "N.n": 0,
                    }
                }
            }
        """
        if not self.proteinGroupsDB:
            return
        self.formattedProteinGroups = {}
        groups: List[ProteinGroup]
        for tableNum, groups in self.proteinGroupsDB.items():
            for group in groups:
                reprAccession = group.representativeAccession
                if reprAccession is None:
                    raise RepresentativeAccessionNotFoundError(group)
                if reprAccession not in self.formattedProteinGroups:
                    self.formattedProteinGroups[reprAccession] = {}
                for accession in group.accessions:
                    if (accession
                            not in self.formattedProteinGroups[reprAccession]):
                        self.formattedProteinGroups[
                            reprAccession][accession] = [
                                Decimal(0) for tableNum in
                                self.proteinGroupsDB.GetSortedTableNums()]
                    self.formattedProteinGroups[reprAccession][accession][
                        self.proteinGroupsDB.GetSortedTableNums().index(
                            tableNum
                        )
                    ] = group.unused

    def GenerateOutputFiles(self) -> None:
        """Создаёт все выходные файлы"""
        if not path.exists(self.inputParams.outputPath):
            makedirs(self.inputParams.outputPath)
        self._accessionsBunch = (
            self.accessionTables.GenerateAccessionsBunchOverAllTables())

        # поле типа bool обозначает нужно ли добавлять столбцы Description и
        # Sequence Length
        fieldsToFiles: Tuple[Tuple[str, str, bool], ...] = (
            ("Counts", "Pep_counts.txt", False),
            ("SeqlenSumm", "Pep_seq_length_summ.txt", False),
            ("ScNormToFileNormRatio", "Sc_norm.txt", False),
            ("ScSumm", "Sc_summ.txt", False),
            ("PSignalNormToFileNormRatio", "Pep_intensity_norm.txt", False),
            ("PSignalSumm", "Pep_intensity_summ.txt", False),
        )

        for field, filename, isAdditionalColumns in fieldsToFiles:
            self.GenerateTableFileByField(field, filename, isAdditionalColumns)

        self.GenerateSequencesFiles(("sequences.fasta", "sequences.txt"))
        self.GenerateGroupsFile("ProteinGroups.txt")
        self.GenerateCountProteinsInGroupsFile("Proteins in groups.txt")
        self.GenerateSettingsFile("settings.txt")
        self.GenerateJointOutputFile("output.txt")

    def GenerateTableFileByField(
            self,
            fieldName: str,
            filename: str,
            isAdditionalColumns: bool) -> None:
        """Создаёт файл с Accession, разбитым по таблицам, где каждой таблице и
        каждому Accession соответствует значение поля fieldName для данного
        Accession в данной таблице

        Args:
            fieldName: имя поля Accession
            filename: имя выходного файла
            isAdditionalColumns: нужно ли добавлять столбцы Description и
                Sequence Length
        """
        with open(self.GetJoinedOutputFilename(filename), 'w') as outFile:
            outFile.write("Accession")
            if isAdditionalColumns:
                outFile.write("\tDescription\tSequence length")
            outFile.write(
                ("\t{}" * len(self.accessionTables.sortedTableNums)).format(
                    *self.accessionTables.sortedTableNums))
            for accession in sorted(self._accessionsBunch.keys()):
                outFile.write(f"\n{accession}")
                if isAdditionalColumns:
                    outFile.write(f"\t{self.seqDB[accession].desc}"
                                  f"\t{self.seqDB[accession].len}")
                for tableNum in self.accessionTables.sortedTableNums:
                    table = self.accessionTables[tableNum]
                    if accession in table:
                        # Создание локального контекста обязательно, иначе
                        # все округления по какой-то причине начинают работать
                        # неправильно! Например, число 0.60000000001
                        # округляется до 0.600001, а не до 0.6
                        with localcontext() as context:
                            context.prec = max(
                                1,
                                7 + Decimal(
                                    table[accession].__dict__[
                                        fieldName]).adjusted())
                            val = str(+table[accession].__dict__[fieldName])
                            if '.' in val and "e" not in val.lower():
                                val = val.rstrip("0").rstrip(".")
                            outFile.write("\t{}".format(val))
                    else:
                        outFile.write('\t')

    def GenerateSequencesFiles(self, filenames: Tuple[str, str]) -> None:
        """Генерирует списки последовательностей по найденным Accession

        Args:
            filenames: имена выходных файлов - первый элемент для fasta файла,
                второй - для txt файла, который можно будет импортировать в
                Excel
        """
        fastaFilename = filenames[0]
        txtFilename = filenames[1]
        accessionsBunch: Set[Sequence] = set()
        for tableName, table in self.accessionTables.items():
            for accessionName in table:
                accessionsBunch.add(accessionName)
        sortedAccessionsBunch = sorted(accessionsBunch)
        with open(self.GetJoinedOutputFilename(fastaFilename), 'w') as outFile:
            for accessionName in sortedAccessionsBunch:
                currentSeqence = self.seqDB[accessionName]
                outFile.write(f'>{accessionName}')
                if len(currentSeqence.desc) > 0:
                    outFile.write(f' {currentSeqence.desc}')
                outFile.write('\n')
                for i in range(0, currentSeqence.len, 60):
                    outFile.write(
                        currentSeqence.seq[i:min(i + 60, currentSeqence.len)])
                    outFile.write('\n')

        with open(self.GetJoinedOutputFilename(txtFilename), 'w') as outFile:
            outFile.write('ID\tSequence\n')
            for accessionName in sortedAccessionsBunch:
                outFile.write(
                    f'{accessionName}\t{self.seqDB[accessionName].seq}\n')

    def GenerateGroupsFile(self, filename: str) -> None:
        """Создаёт файл с Protein группами

        Args:
            filename: имя выходного файла
        """
        if not (self.formattedProteinGroups and self.proteinGroupsDB):
            return
        with open(self.GetJoinedOutputFilename(filename), 'w') as outFile:
            outFile.write("Representative\tAccession" +
                          ("\t{}" * len(
                              self.proteinGroupsDB.GetSortedTableNums())
                           ).format(
                               *self.proteinGroupsDB.GetSortedTableNums()) +
                          "\n")
            for reprAccession, accessions in sorted(
                    self.formattedProteinGroups.items()
            ):
                if len(accessions) == 1:
                    continue
                outFile.write(
                    f"{reprAccession}\t" +
                    ("\t{}" * len(self.proteinGroupsDB.GetSortedTableNums())
                     ).format(*accessions[reprAccession]) + "\n")

                for accession, accessionTables in sorted(accessions.items()):
                    if accession == reprAccession:
                        continue
                    outFile.write(
                        f"\t{accession}" +
                        ("\t{}" * len(
                            self.proteinGroupsDB.GetSortedTableNums())
                         ).format(*accessionTables) + "\n")

    def GenerateCountProteinsInGroupsFile(self, filename: str) -> None:
        """Создаёт файл со всеми репрезентативными Accession и количеством
        Accession в группе

        Репрезентативный Accession не учитывается при подсчёте. Файл имеет вид:
        Accession    ID in group
        RID1         0
        RID2         10
        RID3         3
        RID4         0

        Args:
            filename: имя выходного файла
        """
        if not self.formattedProteinGroups:
            return
        with open(self.GetJoinedOutputFilename(filename), 'w') as outFile:
            outFile.write("Accession\tProteins in group\n")
            for reprAccession, accessions in sorted(
                self.formattedProteinGroups.items()
            ):
                outFile.write(f"{reprAccession}\t{len(accessions) - 1}\n")

    def GenerateSettingsFile(self, filename: str) -> None:
        """Создаёт файл с информацией о параметрах запуска скрипта

        Args:
            filename: имя выходного файла
        """
        with open(self.GetJoinedOutputFilename(filename),
                  'w') as outFile:
            outFile.write(
                '"ProteinPilot summary analyzer"'
                "\n#Protein filter-"
                # TODO
                "\nGlobal FDR critical value (<% k r or default): default"
                f"\nID exclusion list:" + (
                    ' ' + self.inputParams.blackList[0] if (
                        self.inputParams.blackList is not None
                    ) else "") +
                f"\nProtein group filter (Y/N or collapse): " + (
                    "Y" if self.inputParams.isProteinGroupFilter else "N") +
                f"\nPeptide confidence:" + (
                    ' ' + str(self.inputParams.proteinConfidence)
                    if self.inputParams.isProteinConfidence
                    else '')
                + f"\n#Peptide filter-"
                "\nPeptide confidence:"
                f"{(' ' + str(self.inputParams.confPeptide)).strip()}"
                f"\n#Output filter-"
                f"\nMin groups with ID: "
                + str(self.inputParams.minGroupsWithAccession)
                + f"\nMax missing values per group: "
                + str(self.inputParams.maxGroupAbsence)
            )

    def GenerateJointOutputFile(self, filename: str) -> None:
        """Создаёт общий файл с параметрами Accession

        Args:
            filename: имя выходного файла
        """
        with open(self.GetJoinedOutputFilename(filename), 'w') as outFile:
            outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t" +
                          "counts\tSc_summ\tPep_intensity__summ\tSc_norm\t" +
                          "Pep_intensity__norm\tseq_length")
            for accessionName, accessionTables in sorted(
                    self._accessionsBunch.items()):
                for tableNum in sorted(accessionTables.keys(),
                                       key=lambda x: float(x)):
                    accession = accessionTables[tableNum]
                    outFile.write(
                        ("\n{accession}\t{tableNum}\t" +
                         "{seqlenSumm}\t{counts}\t{scSumm}\t" +
                         "{pSignalSumm}\t{scNorm}\t{pSignalNorm}\t" +
                         "{seqlen}").format(
                             accession=accessionName,
                             tableNum=tableNum,
                             seqlenSumm=accession.SeqlenSumm,
                             counts=accession.Counts,
                             scSumm=accession.ScSumm,
                             pSignalSumm=accession.PSignalSumm,
                             scNorm=accession.ScNormToFileNormRatio,
                             pSignalNorm=accession.PSignalNormToFileNormRatio,
                             seqlen=self.seqDB[accessionName].len))

    def GetJoinedOutputFilename(self, filename: str):
        return path.join(self.inputParams.outputPath, filename)
