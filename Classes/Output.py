from typing import Dict, Tuple, List
from os import mkdir, path
from Classes.AccessionTables import AccessionTables
from Classes.ProteinDB import ProteinDB
from Classes.Sequence import Sequence
from Classes.Accession import Accession
from Classes.ProteinGroup import ProteinGroup
from Classes.Errors import AccessionNotFoundError
from decimal import localcontext


class Output:

    outputDirPath: str
    accessionsBunch: Dict[str, Dict[str, Accession]]
    seqDB: Dict[str, Sequence]
    accessionTables: AccessionTables

    def __init__(
            self,
            outputDirPath: str,
            seqDB: Dict[str, Sequence],
            accessionTables: AccessionTables,
            proteinTables: ProteinDB = None) -> None:

        self.seqDB: Dict[str, Sequence] = seqDB
        self.outputDirPath: str = outputDirPath
        self.accessionTables: AccessionTables = accessionTables
        self.GenerateOutputFiles(proteinTables)

    def GenerateOutputFiles(
            self,
            proteinDB: ProteinDB = None) -> None:

        self.CreateDirIfNotExist()
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

        if proteinDB is not None:
            self.GenerateGroupsFile(
                    "ProteinGroups.txt",
                    self.ConvertProteinGroupsToOutputFormat(
                        proteinDB.GetSortedTableNums(),
                        proteinDB.proteinGroupsPerTable),
                    proteinDB)
            self.GenerateGroupsFile(
                    "DifficultCases.txt",
                    self.ConvertProteinGroupsToOutputFormat(
                        proteinDB.GetSortedTableNums(),
                        proteinDB.difficultCases),
                    proteinDB)

        self.GenerateJointOutputFile("output.txt")

    def CreateDirIfNotExist(self) -> None:

        pathLeftover = self.outputDirPath
        folderPath = ""
        while(not path.exists(self.outputDirPath)):
            folderPath += pathLeftover.split('/')[0] + '/'
            pathLeftover = '/'.join(pathLeftover.split('/')[1:])
            try:
                mkdir(folderPath)
            except FileExistsError:
                pass

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
                            context.prec = 6
                            val = str(+table[accession].__dict__[fieldName])
                            if '.' in val and "e" not in val.lower():
                                val = val.rstrip("0").rstrip(".")
                            outFile.write("\t{}".format(val))
                    else:
                        outFile.write('\t')

    def GenerateGroupsFile(
            self,
            outFilename: str,
            groups: Dict[str, Dict[str, Dict[str, int]]],
            proteinDB: ProteinDB) -> None:
        with open(path.join(self.outputDirPath, outFilename), 'w') as outFile:
            outFile.write(("Representative\tAccession" +
                           "\t{}" * len(proteinDB.GetSortedTableNums()) +
                           '\n').format(*proteinDB.GetSortedTableNums()))
            for representativeAccessionName, accessions in sorted(
                    groups.items()):
                if len(accessions) > 1:
                    outFile.write(
                        f"{representativeAccessionName}\t" +
                        ("\t{}" * len(proteinDB.GetSortedTableNums()) +
                         '\n').format(
                             *(map(
                                 lambda key:
                                 accessions[representativeAccessionName][key],
                                 proteinDB.GetSortedTableNums()))
                             ))
                    for replaceableName in sorted(accessions.keys()):
                        if replaceableName == representativeAccessionName:
                            continue
                        outFile.write(
                            ("\t{}" +
                             "\t{}" * len(proteinDB.GetSortedTableNums()) +
                             '\n').format(
                                 replaceableName,
                                 *[accessions[replaceableName][key] for key in
                                     proteinDB.GetSortedTableNums()]
                                 ))

    def ConvertProteinGroupsToOutputFormat(
            self,
            tableNums: List[str],
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
                                    tableNums}
                        )
                    accessions[reprAccession][accession][tableNum] = 1
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
