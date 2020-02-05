from typing import List, Dict, Union, Tuple
from os import listdir, mkdir
from os.path import exists
from Classes import Sequence, Comparable, Accession, Input
from sys import argv

INPUTPATH = "Input"

""" peptideTables - словарь вида
    {
        "1.1": {
            "Заголовок 1": ["Значение 1", "Значение 2", ..., "Значение n1"],
            "Заголовок 2": ["Значение 1", "Значение 2", ..., "Значение n1"],
            ..............................................................
            "Заголовок N1": ["Значение 1", "Значение 2", ..., "Значение n1"]
        },
        "1.2": {
            "Заголовок 1": ["Значение 1", "Значение 2", ..., "Значение n2"],
            "Заголовок 2": ["Значение 1", "Значение 2", ..., "Значение n2"],
            ..............................................................
            "Заголовок N2": ["Значение 1", "Значение 2", ..., "Значение n2"]
        },
        ...,
        "I-ый файл": {
            "Заголовок 1": ["Значение 1", "Значение 2", ..., "Значение nI"],
            "Заголовок 2": ["Значение 1", "Значение 2", ..., "Значение nI"],
            ..............................................................
            "Заголовок NI": ["Значение 1", "Значение 2", ..., "Значение nI"]
        }
    }
"""


def RemoveRow(table: Dict[str, List[str]], rowNum: int) -> None:
    """ Удаляет строку из таблицы

    Удаляет строку из таблицы вида
    {
        "Заголовок 1": ["Значение 1", "Значение 2", ..., "Значение n1"],
        "Заголовок 2": ["Значение 1", "Значение 2", ..., "Значение n1"],
        ..............................................................
        "Заголовок N1": ["Значение 1", "Значение 2", ..., "Значение n1"]
    } """

    columns = [column for column in table]
    for column in columns:
        del table[column][rowNum]


def GenerateJointOutputFile(
        accessionsBunch: Dict[str, Dict[str, Accession]],
        seqDB: Dict[str, Sequence],
        outFilename: str = "output.txt") -> None:

    with open(outFilename, 'w') as outFile:
        outFile.write("Accession\tFilename\tUnused\tseq_length_summ\t\
counts\tSc_summ\tPsignal_summ\tSc_norm\tPsignal_norm\tSP_2\tseq_length")
        for accessionName, accessionTables in accessionsBunch.items():
            for tableNum, accession in accessionTables.items():
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


def GenerateTableFileByField(
        fieldName: str,
        accessionsBunch: Dict[str, Dict[str, Accession]],
        accessionsPerTable: Dict[str, Dict[str, Accession]],
        outFilename: str) -> None:

    with open(outFilename, mode='w') as outFile:
        outFile.write("Accession")
        sortedAccessionsBunchKeys = [*accessionsPerTable]
        sortedAccessionsBunchKeys.sort()
        sortedTableNums = [*accessionsPerTable.keys()]
        sortedTableNums.sort()
        outFile.write((("\t{}" * len(accessionsPerTable))).format(
            *sortedAccessionsBunchKeys))
        for accession in sortedAccessionsBunchKeys:
            outFile.write("\n" + accession)
            for tableNum in sortedTableNums:
                table = accessionsPerTable[tableNum]
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
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> None:

    CreateDirIfNotExist(outputDirPath)
    accessionsBunch = GenerateAccessionsBunchOverAllTables(
        accessionsPerTable)
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
        GenerateTableFileByField(fieldName=field,
                                 accessionsBunch=accessionsBunch,
                                 accessionsPerTable=accessionsPerTable,
                                 outFilename=outputDirPath + '/' + filename)
    GenerateJointOutputFile(accessionsBunch,
                            seqDB,
                            outputDirPath + '/' + "output.txt")


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


def GetAccessionsFromTable(
        peptideTable: Dict[str, List[str]]) -> Dict[str, Accession]:
    """ Получаем unused, суммы значений Sc, Precursor Signal и сумму длинн
    последовательностей и подсчитываем количество строк с одинаковым
    Accession для каждого Accession """

    accessions: Dict[str, Accession] = {}
    i = 0
    while i < len(peptideTable["Accessions"]):
        curAccession = peptideTable["Accessions"][i].split(sep=';')[0]
        if curAccession not in accessions:
            accessions[curAccession] = Accession(name=curAccession)
        accessions[curAccession].Counts += 1
        accessions[curAccession].Unused = float(peptideTable["Unused"][i])
        accessions[curAccession].ScSumm += float(peptideTable["Sc"][i])
        accessions[curAccession].PSignalSumm += float(
            peptideTable["PrecursorSignal"][i])
        accessions[curAccession].SeqlenSumm += len(peptideTable["Sequence"][i])
        i += 1
    return accessions


def CalculateNormParamsForAccessions(accessions: Dict[str, Accession],
                                     seqences: Dict[str, Sequence]):
    for accession in accessions:
        curAccession = accessions[accession]
        curAccession.ScNorm = curAccession.ScSumm / seqences[accession].len
        curAccession.PSignalNorm = (
            curAccession.PSignalSumm / seqences[accession].len)


def GetAccessionsPerTable(
        peptideTables: Dict[str, Dict[str, List[str]]],
        seqences: Dict[str, Sequence]) -> Dict[str, Dict[str, Accession]]:
    """ Получаем суммы значений Sc, Precursor Signal и сумму длинн
    последовательностей для каждого Accession, а также нормализованные
    значения Precursor Signal и Sc для каждого файла"""

    accessionsPerTable: Dict[str, Dict[str, Accession]] = {}
    for tableNum in peptideTables:
        accessionsPerTable[tableNum] = (
            GetAccessionsFromTable(peptideTables[tableNum]))
        CalculateNormParamsForAccessions(accessionsPerTable[tableNum],
                                         seqences)
    return accessionsPerTable


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


def TestConfDefaultCondition(confVal: int):
    if float(confVal) >= 99.0:
        return 2
    if float(confVal) >= 95.0:
        return 1
    return 0


def GenerateTableAccessionsBunch(peptideTable: Dict[str, List[str]]):
    accessionBunch: Dict[str, int] = {}
    for accession in peptideTable["Accessions"]:
        if accession not in accessionBunch:
            accessionBunch[accession.split(';')[0]] = 0
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
            curAccession = curTable["Accessions"][i].split(';')[0]
            if curAccession in blackList:
                blackList[curAccession] += (
                    TestConfDefaultCondition(int(curTable["Conf"][i])))
                if blackList[curAccession] > 1:
                    del blackList[curAccession]
        RemoveAccessionsFromTableByBlacklist(curTable, blackList)
    return peptideTables


def GetFileLines(filename: str) -> Union[List[str], None]:
    """ Returns list of file strings without newline symbols """

    try:
        with open(filename) as tfile:
            return tfile.read().split('\n')
    except FileNotFoundError:
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
    if len(argv) == 9:
        inputParams.whiteList = GetFileLines(argv[5])
        inputParams.blackList = GetFileLines(argv[6])
        inputParams.seqDB = ReadSeqDB(argv[1])
        inputParams.unused = Comparable.GetComparableClass(argv[2])
        inputParams.contrib = Comparable.GetComparableClass(argv[3])
        inputParams.conf = argv[4]
        inputParams.minGroupsWithAccession = int(argv[7])
        inputParams.maxGroupAbsence = int(argv[8])
    else:
        inputParams.whiteList = GetFileLines(input("Id file name: "))
        inputParams.blackList = GetFileLines(input("ID exclusion file name: "))
        inputParams.seqDB = ReadSeqDB(input("Database file name: "))
        inputParams.unused = Comparable.GetComparableClass(input("Unused: "))
        inputParams.contrib = Comparable.GetComparableClass(
            input("Contribution: "))
        inputParams.conf = input("Confidence: ")
        inputParams.minGroupsWithAccession = int(input("Min groups with ID: "))
        inputParams.maxGroupAbsence = int(
            input("Max missing values per group: "))
    return inputParams


def RemoveReversedAccessionsFromTables(
        peptideTables: Dict[str, Dict[str, List[str]]]):

    for peptideTable in peptideTables.values():
        i = 0
        while i < len(peptideTable["Accessions"]):
            if peptideTable["Accessions"][i].startswith("RRRRR"):
                break
            i += 1

        if i < len(peptideTable["Accessions"]):
            for column in peptideTable.values():
                del column[i:]


def ApplyProteinReplacements(replacementsPerTable: Dict[str, Dict[str, str]],
                             peptideTables: Dict[str, Dict[str, List[str]]]):

    for tableName, table in peptideTables.items():
        tableReplacements = replacementsPerTable[tableName]
        for i in range(0, len(table["Accessions"])):
            if table["Accessions"][i] in tableReplacements:
                table["Accessions"][i] = (
                    tableReplacements[table["Accessions"][i]])


def GetRepresentativeAccessionForGroup(
        accessionsWithMaxUnused: Dict[str, float],
        group: List[Tuple[str, float]]) -> str:

    representativeAccession: Tuple[str, float] = (
        group[0][0], accessionsWithMaxUnused[group[0][0]])
    for accessionName, unused in group:
        if accessionsWithMaxUnused[accessionName] > representativeAccession[1]:
            representativeAccession = (accessionName,
                                       accessionsWithMaxUnused[accessionName])
    return representativeAccession[0]


def GetReplacementsForGroups(
        accessionsWithMaxUnused: Dict[str, float],
        groups: List[List[Tuple[str, float]]]
) -> Dict[str, str]:

    replacements = {}
    for group in groups:
        representativeAccession = GetRepresentativeAccessionForGroup(
            accessionsWithMaxUnused, group)
        for accession, unused in group:
            replacements[accession] = representativeAccession
    return replacements


def GetReplacementsPerTable(
        accessionsWithMaxUnused: Dict[str, float],
        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]]
) -> Dict[str, Dict[str, str]]:
    replacementsPerTable: Dict[str, Dict[str, str]] = {}
    for tableName, groups in groupsPerTables.items():
        replacementsPerTable[tableName] = GetReplacementsForGroups(
            accessionsWithMaxUnused, groups)
    return replacementsPerTable


def GetAccessionsWithMaxUnusedFromProteinGroups(
        groupsPerTables: Dict[str, List[List[Tuple[str, float]]]]
) -> Dict[str, float]:
    accessions: Dict[str, float] = {}
    for groups in groupsPerTables.values():
        for group in groups:
            for accessionName, unused in group:
                if accessionName not in accessions:
                    accessions[accessionName] = unused
                elif unused > accessions[accessionName]:
                    accessions[accessionName] = unused
    return accessions


def GetProteinGroupsFromFile(filename: str) -> List[List[Tuple[str, float]]]:
    fileTable = ReadTable(filename)
    groups: List[List[Tuple[str, float]]] = []
    for i in range(0, len(fileTable["N"])):
        if fileTable["Accession"][i].startswith("RRRRR"):
            break
        if float(fileTable["Unused"][i]) != 0:
            unused = float(fileTable["Unused"][i])
            groups.append([])
        groups[-1].append((fileTable["Accession"][i],
                           unused))
    return groups


def GetProteinGroupsFromFiles(
        inputDir: str) -> Dict[str, List[List[Tuple[str, float]]]]:

    groups: Dict[str, List[List[Tuple[str, float]]]] = {}
    for filename in listdir(inputDir):
        if filename.endswith("ProteinSummary.txt"):
            groups[filename.split('_')[0]] = GetProteinGroupsFromFile(
                inputDir + '/' + filename)
    return groups


def GetProteinSummaryReplacements(inputDir: str) -> Dict[str, Dict[str, str]]:
    groupsPerTables = GetProteinGroupsFromFiles(inputDir)
    accessions = GetAccessionsWithMaxUnusedFromProteinGroups(groupsPerTables)
    replacementsPerTable = GetReplacementsPerTable(accessions, groupsPerTables)
    return replacementsPerTable


def ReadTable(tableFilename: str, sep='\t') -> Dict[str, List[str]]:
    """ Считывание файла-таблицы в словарь

    На вход принимается имя файла, из которого считывать таблицу.
    Формат таблицы:
        Заголовок 1\tЗаголовок 2\t...\tЗаголовок N
        Значение 1.1\tЗначение 2.1\t...\tЗначение N.1
        Значение 1.2\tЗначение 2.2\t...\tЗначение N.2
        .............................................
        Значение 1.n\tЗначение 2.n\t...\tЗначение N.n
    Первая строка считается как заголовки столбцов.
    На выходе получается словарь вида:
    {
        "column name 1": ["value1", "value2", ..., "valueN"],
        "column name 2": ["value1", "value2", ..., "valueN"],
        "column name n": ["value1", "value2", ..., "valueN"]
    }
    """

    with open(tableFilename) as tempFile:
        strings = tempFile.read().split('\n')
        tempFile.close()
        table: Dict[str, List[str]] = {}
        columns = strings[0].split(sep)

        for column in columns:
            table[column] = []
        for string in strings[1:]:
            if len(string.strip()):
                i = 0
                for value in string.split(sep):
                    table[columns[i]].append(value)
                    i += 1

        return table
    return None


def ReadPeptideSummaries(inputDir: str) -> Dict:
    """ Считывание всех PeptideSummary файлов в словарь """
    peptideTables = {}
    for filename in listdir(inputDir):
        if "Peptide" in filename:
            peptideTables[filename.split('_')[0]] = (
                ReadTable(inputDir + '/' + filename))
    return peptideTables


def main():
    peptideTables = ReadPeptideSummaries(INPUTPATH)
    proteinReplacements = GetProteinSummaryReplacements(INPUTPATH)
    RemoveReversedAccessionsFromTables(peptideTables)
    ApplyProteinReplacements(proteinReplacements, peptideTables)
    inputParams = GetInput()

    if inputParams.isDefaultConf:
        ApplyDefaultConf(peptideTables)
    ApplyParamsFilter(inputParams.unused,
                      inputParams.contrib,
                      inputParams.conf,
                      peptideTables)

    if inputParams.blackList:
        ApplyBlackList(peptideTables, inputParams.blackList)
    if inputParams.whiteList:
        ApplyWhiteList(peptideTables, inputParams.whiteList)

    accessionsPerFile = GetAccessionsPerTable(peptideTables, inputParams.seqDB)
    filesSumms = GetScPsigAndNormFilesSumm(accessionsPerFile)
    CalculateAccessionsNormRatios(accessionsPerFile, filesSumms)

    ApplyGroupFilter(accessionsPerFile,
                     inputParams.maxGroupAbsence,
                     inputParams.minGroupsWithAccession)

    GenerateOutputFiles("Output/",
                        seqDB=inputParams.seqDB,
                        filesSumms=filesSumms,
                        accessionsPerTable=accessionsPerFile)


if __name__ == "__main__":
    main()
