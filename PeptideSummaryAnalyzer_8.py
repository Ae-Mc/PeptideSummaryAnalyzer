from typing import List, Dict, Union
from os import listdir, mkdir
from os.path import exists
from Classes import Sequence, Comparable, Accession

INPUTPATH = "./Input"

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


def IsFloat(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


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


# TODO: Доделать вывод таблицы в файл
def TableToFile(table: List[List[str]], filename: str):
    with open(filename, mode='w') as outFile:
        for row in table:
            outFile.write('\t'.join(row))


def GenerateTableByField(
        outputDirPath: str,
        fieldName: str,
        accessionsList: List[str],
        accessionsPerTable: Dict[str, Dict[str, Accession]],
        outFilename=None) -> None:

    if outFilename is not None:
        outputFilePath = outputDirPath + '/' + outFilename
    else:
        outputFilePath = outputDirPath + '/' + fieldName.lower() + ".txt"

    with open(outputFilePath, mode='w') as outFile:
        outFile.write("Accession")
        outFile.write((("\t{}" * len(accessionsPerTable)) + '\n').format(
            *accessionsPerTable))
        for i in range(0, len(accessionsList)):
            curAccession = accessionsList[i]
            outFile.write(curAccession)
            for table in accessionsPerTable.values():
                if curAccession in table:
                    outFile.write('\t{}'.format(
                        table[curAccession].__dict__[fieldName]))
                else:
                    outFile.write('\t')
            outFile.write('\n')


def GenerateDescriptionFile(outputDirPath: str,
                            accessionsList: List[str],
                            seqDB: Dict[str, Sequence]) -> None:
    with open(outputDirPath + '/' + "description.txt", mode='w') as descFile:
        descFile.write("Accession\tDescription")
        for accession in accessionsList:
            if len(seqDB[accession].desc):
                descFile.write("\n{}\t{}".format(accession,
                                                 seqDB[accession].desc))


def GenerateAccessionsListOverAllTables(
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> List[str]:

    accessionsList: List[str] = []
    for _, table in accessionsPerTable.items():
        for accession in table:
            if accession not in accessionsList:
                accessionsList.append(accession)
    return accessionsList


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


# TODO: Сделать вывод в файлы
def GenerateOutputFiles(
        outputDirPath: str,
        filesSumms: Dict[str, Dict[str, float]],
        seqDB: Dict[str, Sequence],
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> None:

    CreateDirIfNotExist(outputDirPath)
    accessionsList = GenerateAccessionsListOverAllTables(
        accessionsPerTable)
    GenerateDescriptionFile(outputDirPath, accessionsList, seqDB)
    GenerateTableByField(outputDirPath,
                         fieldName="Counts",
                         accessionsList=accessionsList,
                         accessionsPerTable=accessionsPerTable)
    GenerateTableByField(outputDirPath,
                         fieldName="ScNormToFileNormRatio",
                         accessionsList=accessionsList,
                         accessionsPerTable=accessionsPerTable,
                         outFilename="Sc_norm.txt")
    GenerateTableByField(outputDirPath,
                         fieldName="ScSumm",
                         accessionsList=accessionsList,
                         accessionsPerTable=accessionsPerTable,
                         outFilename="Sc_summ.txt")
    pass


def RemoveAccessionFromAllTables(
        accession: str,
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> None:

    for table in accessionsPerTable.values():
        table.pop(accession, None)


def CountAccessionLackInGroup(
        accession: str,
        group: List[str],
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> int:
    groupAbsence = 0
    for tableName in group:
        if accession not in accessionsPerTable[tableName]:
            groupAbsence += 1
    return groupAbsence


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

    accessionList = GenerateAccessionsListOverAllTables(accessionsPerTable)
    for accession in accessionList:
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
        accessions[curAccession].unused = float(peptideTable["Unused"][i])
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


def GetFileLines(filename: str) -> list:
    """ Returns list of file strings without newline symbols """

    try:
        with open(filename) as tfile:
            return tfile.read().split('\n')
    except FileNotFoundError:
        return None


def GetParam(paramString: str) -> Comparable:
    """ Получение param

    Получение param, общего для всех фалов, либо для каждого своего, из строки.
    Формат: [операция][[имя файла] или [число]]
    Примеры:
        >=99
        < paramList.txt"""

    paramString = paramString.strip()
    param = Comparable(
        op=''.join([ch for ch in paramString if ch in "!=<>"]))

    paramString = paramString[len(param.op):].strip()
    if IsFloat(paramString):
        param.val = float(paramString)
    elif len(param.op) and len(paramString):
        with open(paramString) as paramStringFile:
            strings = paramStringFile.read().replace(' ', '\t').split('\n')
            paramStrings = {}
            for string in strings:
                string = string.strip()
                if len(string):
                    paramStrings[string.split('\t')[0]] = string.split('\t')[1]
            param.val = paramStrings
    return param


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
                    seqDB[seqID].len = len(seqDB[seqID].seq)
                    if not seqDB[seqID].len:
                        input("""Error! Length of sequence with id {} = 0
""".format(seqID))
                        raise(IndexError)
            i += 1

        return seqDB
    return None


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
    seqDB = ReadSeqDB("EFRA_cont.fasta")
    unused = GetParam(">=20")
    contrib = GetParam(">=1")
    conf = ">=50"
    whiteList = GetFileLines("")
    blackList = GetFileLines("IDexcl.txt")

    if conf.strip() == "default":
        ApplyDefaultConf(peptideTables)
        ApplyParamsFilter(unused, contrib, GetParam(""), peptideTables)
    else:
        ApplyParamsFilter(unused, contrib, GetParam(conf), peptideTables)

    minGroupsWithAccession = 1
    maxGroupAbsence = 5

    if blackList:
        ApplyBlackList(peptideTables, blackList)
    if whiteList:
        ApplyWhiteList(peptideTables, whiteList)

    accessionsPerFile = GetAccessionsPerTable(peptideTables, seqDB)
    filesSumms = GetScPsigAndNormFilesSumm(accessionsPerFile)
    CalculateAccessionsNormRatios(accessionsPerFile, filesSumms)

    ApplyGroupFilter(accessionsPerFile,
                     maxGroupAbsence,
                     minGroupsWithAccession)

    GenerateOutputFiles("Output/", seqDB=seqDB, filesSumms=filesSumms,
                        accessionsPerTable=accessionsPerFile)

    for filename in accessionsPerFile:
        continue
        print("{}:".format(filename))
        curFile = accessionsPerFile[filename]
        for accession in curFile:
            curAccession = curFile[accession]
            print("""    name: {name}\n    counts: {counts}
    unused: {unused}\n    ScSumm: {ScSumm}\n    ScNorm: {ScNorm}
    ScNormToFileNormRatio: {ScNormToFile}
    Precursor Signal Summ: {PSignalSumm}
    Precursor Signal Norm: {PSignalNorm}
    PSignalNormToFileNormRatio: {PSignalNormToFile}
    SP/2: {SP_2}
    SeqlenSumm: {SeqlenSumm}\n""".format(
                      name=curAccession.name, counts=curAccession.Counts,
                      unused=curAccession.unused, ScSumm=curAccession.ScSumm,
                      ScNorm=curAccession.ScNorm,
                      ScNormToFile=curAccession.ScNormToFileNormRatio,
                      PSignalSumm=curAccession.PSignalSumm,
                      PSignalNorm=curAccession.PSignalNorm,
                      PSignalNormToFile=(
                          curAccession.PSignalNormToFileNormRatio),
                      SP_2=curAccession.PSignalAndScNormRatiosAverage,
                      SeqlenSumm=curAccession.SeqlenSumm))


if __name__ == "__main__":
    main()
