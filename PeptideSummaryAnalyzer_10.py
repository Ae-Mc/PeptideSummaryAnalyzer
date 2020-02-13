#!python3.7
from typing import List, Dict, Union, Tuple
from os import mkdir
from os.path import exists
from sys import argv
from Classes.Sequence import Sequence
from Classes.Comparable import Comparable
from Classes.Accession import Accession
from Classes.Input import Input
from Classes.ProteinTables import ProteinTables
from Classes.PeptideTables import PeptideTables
from Classes.AccessionTables import AccessionTables

INPUTPATH = "Input"


def RemoveRow(table: Dict[str, List[str]], rowNum: int) -> None:
    columns = [column for column in table]
    for column in columns:
        del table[column][rowNum]


def GenerateJointOutputFile(
        outFilename: str,
        accessionsBunch: Dict[str, Dict[str, Accession]],
        seqDB: Dict[str, Sequence]) -> None:

    with open(outFilename, 'w') as outFile:
        outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t\
counts\tSc_summ\tPsignal_summ\tSc_norm\tPsignal_norm\tSP_2\tseq_length")
        for accessionName, accessionTables in accessionsBunch.items():
            for tableNum in sorted(accessionTables.keys(),
                                   key=lambda x: float(x)):
                accession = accessionTables[tableNum]
                outFile.write("\n{accession}\t{tableNum}\t{unused}\t\
{seqlenSumm}\t{counts}\t{scSumm}\t{pSignalSumm}\t{scNorm}\t{pSignalNorm}\t\
{sp2}\t{seqlen}".format(accession=accessionName,
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


def GenerateGroupsFile(outFilename: str,
                       proteinTables: ProteinTables) -> None:
    with open(outFilename, 'w') as outFile:
        outFile.write(("Representative\tAccession" +
                       "\t{}" * len(proteinTables.sortedTableNums) +
                       '\n').format(*proteinTables.sortedTableNums))
        for representativeAccessionName in \
                proteinTables.proteinReplacementsGroups:
            accession = (
                proteinTables.proteinReplacementsGroups[
                    representativeAccessionName])
            outFile.write(representativeAccessionName)
            for replaceableName in accession:
                outFile.write(
                    ("\t{}" +
                     "\t{}" * len(proteinTables.sortedTableNums) +
                     '\n').format(
                         replaceableName,
                         *[accession[replaceableName][key] for key in
                           proteinTables.sortedTableNums]))


def GenerateTableFileByField(
        fieldName: str,
        accessionsBunch: Dict[str, Dict[str, Accession]],
        accessionTables: AccessionTables,
        outFilename: str) -> None:

    with open(outFilename, mode='w') as outFile:
        outFile.write("Accession")
        outFile.write((("\t{}" * len(accessionTables.sortedTableNums))).format(
            *accessionTables.sortedTableNums))
        for accession in accessionsBunch.keys():
            outFile.write("\n" + accession)
            for tableNum in accessionTables.sortedTableNums:
                table = accessionTables.accessionsPerTable[tableNum]
                if accession in table:
                    outFile.write('\t{}'.format(
                        table[accession].__dict__[fieldName]))
                else:
                    outFile.write('\t')


def GenerateDescriptionFile(outputDirPath: str,
                            accessionsBunch: Dict[str, Dict[str, Accession]],
                            seqDB: Dict[str, Sequence]) -> None:
    with open(outputDirPath + '/' + "description.txt", mode='w') as descFile:
        descFile.write("Accession\tDescription")
        for accession in accessionsBunch.keys():
            if seqDB[accession].len:
                descFile.write("\n{}\t{}".format(accession,
                                                 seqDB[accession].desc))


def GenerateAccessionsBunchOverAllTables(
        accessionsPerTable: Dict[str, Dict[str, Accession]]
) -> Dict[str, Dict[str, Accession]]:

    accessions: Dict[str, Dict[str, Accession]] = {}
    for tableName, table in accessionsPerTable.items():
        for accessionName, accession in table.items():
            if accessionName not in accessions:
                accessions[accessionName] = {}
            accessions[accessionName][tableName] = accession
    return accessions


def CreateDirIfNotExist(directoryPath: str) -> None:

    pathLeftover = directoryPath
    folderPath = ""
    while(not exists(directoryPath)):
        folderPath += pathLeftover.split('/')[0] + '/'
        pathLeftover = '/'.join(pathLeftover.split('/')[1:])
        try:
            mkdir(folderPath)
        except FileExistsError:
            pass


def GenerateOutputFiles(
        outputDirPath: str,
        filesSumms: Dict[str, Dict[str, float]],
        seqDB: Dict[str, Sequence],
        accessionTables: AccessionTables,
        proteinTables: ProteinTables) -> None:

    CreateDirIfNotExist(outputDirPath)
    accessionsBunch = GenerateAccessionsBunchOverAllTables(
        accessionTables.accessionsPerTable)
    GenerateDescriptionFile(outputDirPath, accessionsBunch, seqDB)
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
        GenerateTableFileByField(
            fieldName=field,
            accessionsBunch=accessionsBunch,
            accessionTables=accessionTables,
            outFilename=outputDirPath + '/' + filename)

    GenerateGroupsFile(outputDirPath + '/' + "Groups.txt", proteinTables)

    GenerateJointOutputFile(outputDirPath + '/' + "output.txt",
                            accessionsBunch,
                            # accessionTables=accessionTables,
                            seqDB=seqDB)


def CountAccessionLackInGroup(
        accession: str,
        group: List[str],
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> int:
    groupAbsence = 0
    for tableName in group:
        if accession not in accessionsPerTable[tableName]:
            groupAbsence += 1
    return groupAbsence


def RemoveAccessionFromAllTables(
        accession: str,
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> None:

    for table in accessionsPerTable.values():
        accessionsPerTable.pop(accession, None)


def CountGroupsWithAccession(
        groups: Dict[str, List[str]],
        accession: str,
        maxGroupAbsence: int,
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> int:

    groupsWithAccessionCount = len(groups)
    for tableNames in groups.values():
        if CountAccessionLackInGroup(accession,
                                     tableNames,
                                     accessionsPerTable) > maxGroupAbsence:
            groupsWithAccessionCount -= 1
    return groupsWithAccessionCount


def GenerateGroupsBunch(
        accessionsPerTable: Dict[str, Dict[str, Accession]]
) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}

    for tableNum in accessionsPerTable.keys():
        groupNum = tableNum.split('.')[0]
        if groupNum not in groups:
            groups[groupNum] = []
        groups[groupNum].append(tableNum)
    return groups


def ApplyGroupFilter(accessionsPerTable: Dict[str, Dict[str, Accession]],
                     maxGroupAbsence: int,
                     minGroupsWithAccession: int):
    """ Применяем фильтер по группам

    Таблицы, начинающиеся с одинакового числа входят в одну группу, например:
        1.1, 1.2, 1.3 входят в одну группу, 2.1, 2.2 входят в другую и т. д.
    Из всех таблиц удаляются все Accession, присутствующие в меньше, чем
    minGroupsWithAccession группах. Accession считается отсутствуюющим в
    группе, если он отсутствует в больше, чем maxGroupAbsence таблицах внутри
    группы.
    """

    groups: Dict[str, List[str]] = GenerateGroupsBunch(accessionsPerTable)

    accessionBunch = GenerateAccessionsBunchOverAllTables(accessionsPerTable)
    for accession in accessionBunch:
        if CountGroupsWithAccession(
                groups,
                accession,
                maxGroupAbsence,
                accessionsPerTable) < minGroupsWithAccession:
            RemoveAccessionFromAllTables(accession, accessionsPerTable)


def CalculateAccessionsNormRatios(
        accessionsPerTable: Dict[str, Dict[str, Accession]],
        tableSumms: Dict[str, Dict[str, float]]) -> None:
    """ NormRatio — отношение Sc/PSignalNorm Accession к сумме всех
    Sc/PsignalNorm в файле """

    for tableNum in tableSumms:
        curSumm = tableSumms[tableNum]
        curAccessionTable = accessionsPerTable[tableNum]
        for accession in curAccessionTable:
            curAccession = curAccessionTable[accession]
            curAccession.ScNormToFileNormRatio = (
                curAccession.ScNorm / curSumm["ScNorm"])
            curAccession.PSignalNormToFileNormRatio = (
                curAccession.PSignalNorm /
                curSumm["PSignalNorm"])
            curAccession.PSignalAndScNormRatiosAverage = (
                (curAccession.ScNormToFileNormRatio +
                 curAccession.PSignalNormToFileNormRatio) / 2)


def GetScPsigAndNormFilesSumm(
    accessionsPerTable: Dict[str, Dict[str, Accession]]) -> Dict[
                        str, Dict[str, float]]:
    """ Получаем суммы параметров Sc, Sequence, PrecursorSignal,
    ScNorm, PSignalNorm по файлам


    Возвращаем словарь с суммами вида
    {
        "1.1": {
            "ScSumm": S1,
            "pSignalSumm": P1,
            "ScNorm: SN1,
            "PSignalNorm": PN1
        },
        "1.2": {
            "ScSumm": S2,
            "pSignalSumm": P2,
            "ScNorm: SN2,
            "PSignalNorm": PN2
        },
        ...,
        "n-ый файл": {
            "ScSumm": Sn,
            "pSignalSumm": Pn,
            "ScNorm: SNn,
            "PSignalNorm": PNn
        }
    }
    """

    fileSumms: Dict[str, Dict[str, Union[float, int]]] = {}

    for tableNum in accessionsPerTable:
        curTable = accessionsPerTable[tableNum]
        fileSumms[tableNum] = {}
        curSumm = fileSumms[tableNum]
        curSumm["ScSumm"] = 0
        curSumm["PSignalSumm"] = 0
        curSumm["ScNorm"] = 0
        curSumm["PSignalNorm"] = 0

        for accession in curTable:
            curSumm["ScSumm"] += float(curTable[accession].ScSumm)
            curSumm["ScNorm"] += float(curTable[accession].ScNorm)
            curSumm["PSignalSumm"] += float(curTable[accession].PSignalSumm)
            curSumm["PSignalNorm"] += float(curTable[accession].PSignalNorm)
    return fileSumms


def ApplyWhiteList(peptideTables: Dict[str, Dict[str, List[str]]],
                   whiteList: List[str]):
    """ Удаляем из всех таблиц все id, отсутствующие в белом списке """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Accessions"])
        i = 0
        while i < curTableLen:
            if curTable["Accessions"][i].split(';')[0] not in whiteList:
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1


def ApplyBlackList(peptideTables: Dict[str, Dict[str, List[str]]],
                   blackList: List[str]):
    """ Удаляем из всех таблиц все id, находящиеся в чёрном списке """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Accessions"])
        i = 0
        while i < curTableLen:
            if curTable["Accessions"][i].split(';')[0] in blackList:
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1


def TestUnusedContribConfParams(unused: Comparable,
                                contrib: Comparable,
                                conf: Comparable,
                                peptideTable: Dict[str, List[str]],
                                tableNum: str,
                                tableRowNum: int):
    if(unused.compare(peptideTable["Unused"][tableRowNum], tableNum) and
       contrib.compare(peptideTable["Contrib"][tableRowNum], tableNum) and
       conf.compare(peptideTable["Conf"][tableRowNum], tableNum)):
        return True
    return False


def ApplyParamsFilter(unused: Comparable,
                      contrib: Comparable,
                      conf: Comparable,
                      peptideTables: Dict[str, Dict[str, List[str]]]) -> dict:
    """ Применяем фильтры, завязанные на параметры unused, contib, conf.

    Возвращаем словарь с уже применёнными фильрами.
    """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Unused"])
        i = 0

        while i < curTableLen:
            if not TestUnusedContribConfParams(unused, contrib, conf,
                                               curTable, tableNum, i):
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1

    return peptideTables


def RemoveAccessionsFromTableByBlacklist(peptideTable: Dict[str, List[str]],
                                         blackList: Dict[str, int]) -> None:
    i = 0
    curTableLen = len(peptideTable["Accessions"])
    while i < curTableLen:
        if(peptideTable["Accessions"][i].split(';')[0] in blackList):
            RemoveRow(peptideTable, i)
            curTableLen -= 1
            i -= 1
        i += 1


def TestConfDefaultCondition(confVal: float):
    if confVal >= 99.0:
        return 2
    if confVal >= 95.0:
        return 1
    return 0


def GenerateTableAccessionsBunch(peptideTable: Dict[str, List[str]]):
    accessionBunch: Dict[str, int] = {}
    for accession in peptideTable["Accessions"]:
        if accession not in accessionBunch:
            accessionBunch[accession] = 0
    return accessionBunch


def ApplyDefaultConf(peptideTables: Dict[str, Dict[str, List[str]]]):
    """ default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95"""

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Unused"])
        # Заносим все Accession в список Accession для удаления
        # В процессе чтения списка записи из него будут удаляться
        blackList = GenerateTableAccessionsBunch(curTable)
        i = 0
        while i < curTableLen:
            curAccession = curTable["Accessions"][i]
            if curAccession in blackList:
                blackList[curAccession] += (
                    TestConfDefaultCondition(float(curTable["Conf"][i])))
                if blackList[curAccession] > 1:
                    del blackList[curAccession]
            i += 1
        RemoveAccessionsFromTableByBlacklist(curTable, blackList)
    return peptideTables


def GetFileLines(filename: str) -> Union[List[str], None]:
    """ Returns list of file strings without newline symbols """

    if len(filename):
        with open(filename) as tfile:
            return tfile.read().split('\n')
    return None


def ReadSeqDB(seqDBFilename: str) -> Dict[str, Sequence]:
    """ Считывание последовательностей из файла

    Считывание длин последовательностей из файла БД с последовательностями в
    словарь классов Sequence вида {"Accession"} = Sequence """

    with open(seqDBFilename) as seqDBFile:
        strings = seqDBFile.read().split('\n')
        seqDBFile.close()
        seqDB = {}
        i = 0
        while(i < len(strings)):
            if len(strings[i]):
                if strings[i][0] == '>':
                    seqID = strings[i].split(' ')[0][1:]
                    seqDB[seqID] = Sequence()
                    if len(strings[i].split(' ')) > 1:
                        seqDB[seqID].desc = strings[i].split(' ')[1]
                        for word in strings[i].split(' ')[2:]:
                            if word.startswith("OS="):
                                break
                            seqDB[seqID].desc += ' ' + word
                    i += 1
                    while ((i < len(strings))
                           and (not len(strings[i]) or strings[i][0] != '>')):
                        seqDB[seqID].seq += ''.join(
                            [ch for ch in strings[i]
                             if ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"])
                        i += 1
                    i -= 1
                    if not seqDB[seqID].len:
                        input("""Error! Length of sequence with id {} = 0
""".format(seqID))
                        raise(IndexError)
            i += 1

        return seqDB
    return None


def GetInput() -> Input:
    inputParams = Input()
    if len(argv) == 10:
        inputParams.whiteList = GetFileLines(argv[1])
        inputParams.blackList = GetFileLines(argv[2])
        inputParams.isProteinGroupFilter = (
            True if argv[3].strip() == 'y' else False)
        inputParams.seqDB = ReadSeqDB(argv[4])
        inputParams.unused = Comparable(argv[5])
        inputParams.contrib = Comparable(argv[6])
        inputParams.conf = argv[7]
        inputParams.minGroupsWithAccession = int(argv[8])
        inputParams.maxGroupAbsence = int(argv[9])
    else:
        inputParams.whiteList = GetFileLines(input("Id file name: "))
        inputParams.blackList = GetFileLines(input("ID exclusion file name: "))
        inputParams.isProteinGroupFilter = True if input(
            "Protein group filter: ").strip() == 'y' else False
        inputParams.seqDB = ReadSeqDB(input("Database file name: "))
        inputParams.unused = Comparable(input("Unused: "))
        inputParams.contrib = Comparable(input("Contribution: "))
        inputParams.conf = input("Confidence: ")
        inputParams.minGroupsWithAccession = int(input("Min groups with ID: "))
        inputParams.maxGroupAbsence = int(
            input("Max missing values per group: "))
    return inputParams


def main():
    peptideTables = PeptideTables(INPUTPATH)
    proteinTables = ProteinTables(inputDir=INPUTPATH)
    peptideTables.ApplyProteinReplacements(proteinTables)
    inputParams = GetInput()

    if inputParams.isDefaultConf:
        ApplyDefaultConf(peptideTables.peptideTables)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.contrib,
                      inputParams.conf,
                      peptideTables.peptideTables)

    if inputParams.whiteList:
        ApplyWhiteList(peptideTables.peptideTables, inputParams.whiteList)

    accessionTables = AccessionTables(inputParams.seqDB, peptideTables)
    accessionTables.sortedTableNums = peptideTables.GetSortedTableNums()
    filesSumms = GetScPsigAndNormFilesSumm(accessionTables.accessionsPerTable)

    if inputParams.blackList:
        ApplyBlackList(peptideTables.peptideTables, inputParams.blackList)

    accessionTables.GetAccessionsPerTable(inputParams.seqDB, peptideTables)
    CalculateAccessionsNormRatios(accessionTables.accessionsPerTable,
                                  filesSumms)

    ApplyGroupFilter(accessionTables.accessionsPerTable,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    GenerateOutputFiles("Output/",
                        seqDB=inputParams.seqDB,
                        filesSumms=filesSumms,
                        accessionTables=accessionTables,
                        proteinTables=proteinTables)


if __name__ == "__main__":
    main()
