from typing import Dict, Tuple, List, Optional
from os import makedirs, path
from Classes.AccessionTables import AccessionTables
from Classes.Sequence import Sequence
from Classes.Accession import Accession
from Classes.ProteinGroupsDB import ProteinGroupsDB
from Classes.ProteinGroup import ProteinGroup
from Classes.Errors import AccessionNotFoundError
from decimal import localcontext, Decimal


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
    outputDirPath: str
    _accessionsBunch: Dict[str, Dict[str, Accession]]
    seqDB: Dict[str, Sequence]
    accessionTables: AccessionTables
    proteinGroupsDB: Optional[ProteinGroupsDB]

    def __init__(
            self,
            outputDirPath: str,
            seqDB: Dict[str, Sequence],
            accessionTables: AccessionTables,
            proteinGroupsDB: ProteinGroupsDB = None) -> None:

        self.seqDB: Dict[str, Sequence] = seqDB
        self.outputDirPath: str = outputDirPath
        self.accessionTables: AccessionTables = accessionTables
        self.proteinGroupsDB = proteinGroupsDB
        self.GenerateOutputFiles()

    def GenerateOutputFiles(self) -> None:
        """Создаёт все выходные файлы"""
        if not path.exists(self.outputDirPath):
            makedirs(self.outputDirPath)
        self._accessionsBunch = (
            self.accessionTables.GenerateAccessionsBunchOverAllTables())
        self.GenerateDescriptionFile(filename="description.txt")

        fieldsToFiles: Tuple[Tuple[str, str, bool], ...] = (
            ("Counts", "counts.txt", False),
            ("ScNormToFileNormRatio", "Sc_norm.txt", True),
            ("ScSumm", "Sc_summ.txt", True),
            ("PSignalNormToFileNormRatio", "Pep_intensity_norm.txt", True),
            ("PSignalSumm", "Pep_intensity_summ.txt", True),
            ("PSignalAndScNormRatiosAverage", "SP_2.txt", False),
            ("SeqlenSumm", "seq_length_summ.txt", False),
            ("Unused", "unused.txt", False)
        )

        for field, filename, isSeqlen in fieldsToFiles:
            self.GenerateTableFileByField(field, filename, isSeqlen)

        self.GenerateGroupsFile("ProteinGroups.txt")
        self.GenerateJointOutputFile("output.txt")

    def GenerateDescriptionFile(self, filename: str) -> None:
        """Создаёт файл с описаниями каждого Accession

        Args:
            filename: имя выходного файла
        """
        with open(path.join(self.outputDirPath, filename), 'w') as descFile:
            descFile.write("Accession\tDescription")
            for accession in sorted(self._accessionsBunch.keys()):
                if self.seqDB[accession].len:
                    descFile.write("\n{}\t{}".format(
                        accession, self.seqDB[accession].desc))

    def GenerateTableFileByField(
            self,
            fieldName: str,
            filename: str,
            isSeqlen: bool) -> None:
        """Создаёт файл с Accession, разбитым по таблицам, где каждой таблице и
        каждому Accession соответствует значение поля fieldName для данного
        Accession в данной таблице

        Args:
            fieldName: имя поля Accession
            filename: имя выходного файла
            isSeqlen: нужно ли добавлять столбец с Sequence Length
        """
        with open(path.join(self.outputDirPath, filename),
                  'w') as outFile:
            outFile.write("Accession")
            if isSeqlen:
                outFile.write("\tSequence length")
            outFile.write(
                ("\t{}" * len(self.accessionTables.sortedTableNums)).format(
                    *self.accessionTables.sortedTableNums))
            for accession in sorted(self._accessionsBunch.keys()):
                outFile.write(f"\n{accession}")
                if isSeqlen:
                    outFile.write(f"\t{self.seqDB[accession].len}")
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

    def GenerateGroupsFile(self, filename: str) -> None:
        """Создаёт файл с Protein группами

        Args:
            filename: имя выходного файла
        """
        if self.proteinGroupsDB is None:
            return
        outDict = self.ConvertProteinGroupsToOutputFormat(self.proteinGroupsDB)
        with open(path.join(self.outputDirPath, filename), 'w') as outFile:
            outFile.write("Representative\tAccession" +
                          ("\t{}" * len(
                              self.proteinGroupsDB.GetSortedTableNums())
                           ).format(
                               *self.proteinGroupsDB.GetSortedTableNums()) +
                          "\n")
            for reprAccession, accessions in sorted(outDict.items()):
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

    def ConvertProteinGroupsToOutputFormat(
            self, proteinGroupsDB: ProteinGroupsDB
    ) -> Dict[str, Dict[str, List[int]]]:
        """Представление групп в виде, удобном для вывода {
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
            }

        Args:
            proteinGroupsDB: база данных Protein групп

        Returns:
            Представление групп в виде, удобном для вывода: {
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
            }
        """
        accessions: Dict[str, Dict[str, List[int]]] = {}
        groups: List[ProteinGroup]
        for tableNum, groups in proteinGroupsDB.items():
            for group in groups:
                reprAccession = group.representativeAccession
                if reprAccession is None:
                    raise AccessionNotFoundError(
                        "Representative accession for group"
                        f"{group.accessions} not found")
                if reprAccession not in accessions:
                    accessions[reprAccession] = {}
                for accession in group.accessions:
                    if (accession not in accessions[reprAccession]):
                        accessions[reprAccession][accession] = (
                                [0 for tableNum in
                                    proteinGroupsDB.GetSortedTableNums()])
                    accessions[reprAccession][accession][
                        proteinGroupsDB.GetSortedTableNums(
                        ).index(tableNum)] = 1
        return accessions

    def GenerateJointOutputFile(self, filename: str) -> None:
        """Создаёт общий файл с параметрами Accession

        Args:
            filename: имя выходного файла
        """
        with open(path.join(self.outputDirPath, filename), 'w') as outFile:
            outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t" +
                          "counts\tSc_summ\tPep_intensity__summ\tSc_norm\t" +
                          "Pep_intensity__norm\tSP_2\tseq_length")
            for accessionName, accessionTables in sorted(
                    self._accessionsBunch.items()):
                for tableNum in sorted(accessionTables.keys(),
                                       key=lambda x: float(x)):
                    accession = accessionTables[tableNum]
                    outFile.write(
                        ("\n{accession}\t{tableNum}\t{unused}\t" +
                         "{seqlenSumm}\t{counts}\t{scSumm}\t" +
                         "{pSignalSumm}\t{scNorm}\t{pSignalNorm}\t" +
                         "{sp2}\t{seqlen}").format(
                             accession=accessionName,
                             tableNum=tableNum,
                             unused=accession.Unused,
                             seqlenSumm=accession.SeqlenSumm,
                             counts=accession.Counts,
                             scSumm=accession.ScSumm,
                             pSignalSumm=accession.PSignalSumm,
                             scNorm=accession.ScNormToFileNormRatio,
                             pSignalNorm=accession.PSignalNormToFileNormRatio,
                             sp2=accession.PSignalAndScNormRatiosAverage,
                             seqlen=self.seqDB[accessionName].len))
