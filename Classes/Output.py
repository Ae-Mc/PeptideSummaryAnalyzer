from typing import Dict, Tuple
from os import mkdir
from os.path import exists
from Classes.AccessionTables import AccessionTables
from Classes.ProteinTables import ProteinTables
from Classes.Sequence import Sequence
from Classes.Accession import Accession


class Output:

    outputDirPath: str
    accessionsBunch: Dict[str, Dict[str, Accession]]

    def __init__(
            self,
            outputDirPath: str = None,
            filesSumms: Dict[str, Dict[str, float]] = None,
            seqDB: Dict[str, Sequence] = None,
            accessionTables: AccessionTables = None,
            proteinTables: ProteinTables = None) -> None:

        if(outputDirPath is not None and
           filesSumms is not None and
           seqDB is not None and
           accessionTables is not None and
           proteinTables is not None):
            self.GenerateOutputFiles(outputDirPath,
                                     filesSumms,
                                     seqDB,
                                     accessionTables,
                                     proteinTables)

    def GenerateOutputFiles(
            self,
            outputDirPath: str,
            filesSumms: Dict[str, Dict[str, float]],
            seqDB: Dict[str, Sequence],
            accessionTables: AccessionTables,
            proteinTables: ProteinTables) -> None:

        self.outputDirPath = outputDirPath

        self.CreateDirIfNotExist()
        self.accessionsBunch = (
            accessionTables.GenerateAccessionsBunchOverAllTables())
        self.GenerateDescriptionFile(seqDB)

        fieldsToFiles: Tuple[Tuple[str, str], ...] = (
            ("Counts", "counts.txt"),
            ("ScNormToFileNormRatio", "Sc_norm.txt"),
            ("ScSumm", "Sc_summ.txt"),
            ("PSignalNormToFileNormRatio", "Psignal_norm.txt"),
            ("PSignalSumm", "Psignal_summ.txt"),
            ("PSignalAndScNormRatiosAverage", "SP_2.txt"),
            ("SeqlenSumm", "seq_length_summ.txt"),
            ("Unused", "unused.txt")
        )
        for field, filename in fieldsToFiles:
            self.GenerateTableFileByField(
                fieldName=field,
                accessionTables=accessionTables,
                outFilename=filename)

        self.GenerateGroupsFile("Groups.txt", proteinTables)

        self.GenerateJointOutputFile("output.txt",
                                     seqDB=seqDB)

    def CreateDirIfNotExist(self) -> None:

        pathLeftover = self.outputDirPath
        folderPath = ""
        while(not exists(self.outputDirPath)):
            folderPath += pathLeftover.split('/')[0] + '/'
            pathLeftover = '/'.join(pathLeftover.split('/')[1:])
            try:
                mkdir(folderPath)
            except FileExistsError:
                pass

    def GenerateDescriptionFile(
            self,
            seqDB: Dict[str, Sequence],
            outFilename: str = "description.txt") -> None:
        with open(self.outputDirPath +
                  '/' +
                  outFilename,
                  mode='w') as descFile:
            descFile.write("Accession\tDescription")
            for accession in sorted(self.accessionsBunch.keys()):
                if seqDB[accession].len:
                    descFile.write("\n{}\t{}".format(accession,
                                                     seqDB[accession].desc))

    def GenerateTableFileByField(
            self,
            fieldName: str,
            accessionTables: AccessionTables,
            outFilename: str) -> None:

        with open(self.outputDirPath + '/' + outFilename, mode='w') as outFile:
            outFile.write("Accession")
            outFile.write((("\t{}" *
                            len(accessionTables.sortedTableNums))).format(
                *accessionTables.sortedTableNums))
            for accession in sorted(self.accessionsBunch.keys()):
                outFile.write("\n" + accession)
                for tableNum in accessionTables.sortedTableNums:
                    table = accessionTables.accessionsPerTable[tableNum]
                    if accession in table:
                        outFile.write('\t{}'.format(
                            table[accession].__dict__[fieldName]))
                    else:
                        outFile.write('\t')

    def GenerateGroupsFile(
            self,
            outFilename: str,
            proteinTables: ProteinTables) -> None:
        with open(self.outputDirPath + '/' + outFilename, 'w') as outFile:
            outFile.write(("Representative\tAccession" +
                           "\t{}" * len(proteinTables.sortedTableNums) +
                           '\n').format(*proteinTables.sortedTableNums))
            for representativeAccessionName in sorted(
                    proteinTables.proteinReplacementsGroups):
                accession: Dict[str, Dict[str, int]] = (
                    proteinTables.proteinReplacementsGroups[
                        representativeAccessionName])
                outFile.write(representativeAccessionName)
                for replaceableName in sorted(accession.keys()):
                    outFile.write(
                        ("\t{}" +
                         "\t{}" * len(proteinTables.sortedTableNums) +
                         '\n').format(
                             replaceableName,
                             *[accession[replaceableName][key] for key in
                               proteinTables.sortedTableNums]))

    def GenerateJointOutputFile(
            self,
            outFilename: str,
            seqDB: Dict[str, Sequence]) -> None:

        with open(self.outputDirPath + '/' + outFilename, 'w') as outFile:
            outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t" +
                          "counts\tSc_summ\tPsignal_summ\tSc_norm\t" +
                          "Psignal_norm\tSP_2\tseq_length")
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
                             seqlen=seqDB[accessionName].len))
