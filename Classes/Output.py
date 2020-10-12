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

    outputDirPath: str
    accessionsBunch: Dict[str, Dict[str, Accession]]
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

        if not path.exists(self.outputDirPath):
            makedirs(self.outputDirPath)
        self.accessionsBunch = (
            self.accessionTables.GenerateAccessionsBunchOverAllTables())
        self.GenerateDescriptionFile(outFilename="description.txt")

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
            if isSeqlen:
                self.GenerateTableFileByFieldWithSeqlen(
                    fieldName=field,
                    outFilename=filename)
            else:
                self.GenerateTableFileByField(
                    fieldName=field,
                    outFilename=filename)

        self.GenerateGroupsFile("ProteinGroups.txt")
        self.GenerateJointOutputFile("output.txt")

    def GenerateDescriptionFile(
            self,
            outFilename: str) -> None:
        with open(path.join(self.outputDirPath, outFilename), 'w') as descFile:
            descFile.write("Accession\tDescription")
            for accession in sorted(self.accessionsBunch.keys()):
                if self.seqDB[accession].len:
                    descFile.write("\n{}\t{}".format(
                        accession, self.seqDB[accession].desc))

    def GenerateTableFileByFieldWithSeqlen(
            self,
            fieldName: str,
            outFilename: str) -> None:
        with open(path.join(self.outputDirPath, outFilename),
                  'w') as outFile:
            outFile.write(
                "Accession\tSequence length" +
                ("\t{}" * len(self.accessionTables.sortedTableNums)).format(
                    *self.accessionTables.sortedTableNums))
            for accession in sorted(self.accessionsBunch.keys()):
                outFile.write(f"\n{accession}" +
                              f"\t{self.seqDB[accession].len}")
                for tableNum in self.accessionTables.sortedTableNums:
                    table = self.accessionTables.accessionsPerTable[tableNum]
                    if accession in table:
                        # Создание локального контекста обязательно, иначе
                        # все округления по какой-то причине начинают работать
                        # неправильно! Например, число 0.60000000001
                        # округляется до 0.600001, а не до 0.6
                        with localcontext() as context:
                            context.prec = max(
                                1,
                                7 + table[accession
                                          ].__dict__[fieldName].adjusted())
                            val = str(+table[accession].__dict__[fieldName])
                            if '.' in val and "e" not in val.lower():
                                val = val.rstrip("0").rstrip(".")
                            outFile.write("\t{}".format(val))
                    else:
                        outFile.write('\t')

    def GenerateTableFileByField(
            self,
            fieldName: str,
            outFilename: str) -> None:

        with open(path.join(self.outputDirPath, outFilename), 'w') as outFile:
            outFile.write("Accession")
            outFile.write(
                ("\t{}" * len(self.accessionTables.sortedTableNums)).format(
                    *self.accessionTables.sortedTableNums))
            for accession in sorted(self.accessionsBunch.keys()):
                outFile.write("\n" + accession)
                for tableNum in self.accessionTables.sortedTableNums:
                    table = self.accessionTables.accessionsPerTable[tableNum]
                    if accession in table:
                        # Создание локального контекста обязательно, иначе
                        # все округления по какой-то причине начинают работать
                        # неправильно! Например, число 0.60000000001
                        # округляется до 0.600001, а не до 0.6
                        with localcontext() as context:
                            context.prec = max(
                                1,
                                7 + Decimal(
                                        table[accession].__dict__[fieldName]
                                ).adjusted())
                            val = str(
                                +Decimal(table[accession].__dict__[fieldName]))
                            if '.' in val and "e" not in val.lower():
                                val = val.rstrip("0").rstrip(".")
                            outFile.write("\t{}".format(val))
                    else:
                        outFile.write('\t')

    def GenerateGroupsFile(self, outFilename: str) -> None:
        if self.proteinGroupsDB is None:
            return
        outDict = self.ConvertProteinGroupsToOutputFormat(self.proteinGroupsDB)
        with open(path.join(self.outputDirPath, outFilename), 'w') as outFile:
            outFile.write("Representative\tAccession" +
                          ("\t{}" * len(self.proteinGroupsDB.sortedTableNums)
                           ).format(*self.proteinGroupsDB.sortedTableNums) +
                          "\n")
            for reprAccession, accessions in sorted(outDict.items()):
                if len(accessions) == 1:
                    continue
                outFile.write(
                    f"{reprAccession}\t" +
                    ("\t{}" * len(self.proteinGroupsDB.sortedTableNums)
                     ).format(*accessions[reprAccession]) + "\n")

                for accession, accessionTables in sorted(accessions.items()):
                    if accession == reprAccession:
                        continue
                    outFile.write(
                        f"\t{accession}" +
                        ("\t{}" * len(self.proteinGroupsDB.sortedTableNums)
                         ).format(*accessionTables) + "\n")

    def ConvertProteinGroupsToOutputFormat(
            self, proteinGroupsDB: ProteinGroupsDB
    ) -> Dict[str, Dict[str, List[int]]]:
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
                                    proteinGroupsDB.sortedTableNums])
                    accessions[reprAccession][accession][
                        proteinGroupsDB.sortedTableNums.index(tableNum)] = 1
        return accessions

    def GenerateJointOutputFile(
            self,
            outFilename: str) -> None:

        with open(path.join(self.outputDirPath, outFilename), 'w') as outFile:
            outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t" +
                          "counts\tSc_summ\tPep_intensity__summ\tSc_norm\t" +
                          "Pep_intensity__norm\tSP_2\tseq_length")
            for accessionName, accessionTables in sorted(
                    self.accessionsBunch.items()):
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
