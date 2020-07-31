from sys import argv
from typing import List, Dict, Union
from Classes.Accession import Accession
from Classes.AccessionTables import AccessionTables
from Classes.ColumnNames import ColumnNames
from Classes.Comparable import Comparable
from Classes.Input import Input
from Classes.Sequence import Sequence
from decimal import Decimal


## Удаляет строку из таблицы, полученной с помощью Classes.ReadTable.ReadTable
# @param table Таблица, из которой происходит удаление
# @param rowNum Номер строки в таблице
def RemoveRow(table: Dict[str, List[str]], rowNum: int) -> None:
    columns = [column for column in table]
    for column in columns:
        del table[column][rowNum]


## Подсчитывает количество файлов, не содержащих данный accession, в группе
# @param accession Имя accession, для которого считается отсутствие в группе
# @param group Список номеров таблиц, входящих в группу
# @param accessionsPerTable Словарь, ключом в котором является номер таблицы, а
#                           ключом — словарь, в котором, в свою очередь, ключом
#                           является имя accession, а значением — экземпляр
#                           класса Classes.Accession
# @returns Количество файлов, не содержащих данный accession, в группе
def CountAccessionLackInGroup(
        accession: str,
        group: List[str],
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> int:
    groupAbsence = 0
    for tableName in group:
        if accession not in accessionsPerTable[tableName]:
            groupAbsence += 1
    return groupAbsence


## Подсчёт количества групп, в которых данный accession отсутствует не более
#  maxGroupAbsence раз
# @param groups Словарь, в котором ключом является номер группы, а значением —
#               список номеров входящих в неё таблиц
# @param accession Имя accession, для которого считается отсутствие в группах
# @param maxGroupAbsence Максимальное количество отсутствий accession в группе
# @param accessionsPerTable Словарь, ключом в котором является номер таблицы, а
#                           ключом — словарь, в котором, в свою очередь, ключом
#                           является имя accession, а значением — экземпляр
#                           класса Classes.Accession
# @returns Количество групп, в которых accession отсутствует не более
#          maxGroupAbsence раз
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


## Создаёт словарь с группами, в котором ключом является номер группы, а
#  значением — список номеров таблиц, входящих в группу
# @param accessionsPerTable Словарь, ключом в котором является номер таблицы, а
#                           ключом — словарь, в котором, в свою очередь, ключом
#                           является имя accession, а значением — экземпляр
#                           класса Classes.Accession
# @returns Словарь с группами, в котором ключом является номер группы, а
#  значением — список номеров таблиц, входящих в группу
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


## @brief Применение фильтра по группам к экземпляру класса
#  Classes.AccessionTables
# @param accessionTables Экземпляр класса Classes.AccessionTables, к которому
#                        применяется фильтр
#
# Таблицы, начинающиеся с одинакового числа входят в одну группу, например:
#     1.1, 1.2, 1.3 входят в одну группу, 2.1, 2.2 входят в другую и т. д.
# Из всех таблиц удаляются все Accession, присутствующие в меньше, чем в
# minGroupsWithAccession группах. Accession считается отсутствуюющим в
# группе, если он отсутствует в больше, чем maxGroupAbsence таблицах внутри
# группы.
def ApplyGroupFilter(accessionTables: AccessionTables,
                     maxGroupAbsence: int,
                     minGroupsWithAccession: int) -> None:
    groups: Dict[str, List[str]] = GenerateGroupsBunch(
        accessionTables.accessionsPerTable)

    accessionBunch = accessionTables.GenerateAccessionsBunchOverAllTables()
    for accession in accessionBunch:
        if CountGroupsWithAccession(
                groups,
                accession,
                maxGroupAbsence,
                accessionTables.accessionsPerTable) < minGroupsWithAccession:
            accessionTables.RemoveAccessionFromAllTables(accession)


## Подсчёт нормализованных значений PSignalNormToFileNormRatio и
#  ScNormToFileNormRatio для всех accession
# @param accessionTables Экземпляр класса Classes.AccessionTables, к которому
#                        применяется фильтр
# @param tableSumms Словарь, содержащий суммы значений ScNorm и PSignalNorm
#                   для каждой таблицы
#
# Нормализованные значения PSignalNormToFileNormRatio и ScNormToFileNormRatio
# это отношения ScNorm/PsignalNorm accession к сумме всех ScNorm/PsignalNorm в
# файле
def CalculateAccessionsNormRatios(
        accessionTables: AccessionTables,
        tableSumms: Dict[str, Dict[str, Decimal]]) -> None:
    for tableNum, curSumm in tableSumms.items():
        curAccessionTable = accessionTables.accessionsPerTable[tableNum]
        for accession, curAccession in curAccessionTable.items():
            curAccession.ScNormToFileNormRatio = (
                (curAccession.ScNorm /  # noqa: W504
                 curSumm["ScNorm"]) if curSumm["ScNorm"] != 0
                else Decimal(0))
            curAccession.PSignalNormToFileNormRatio = (
                (curAccession.PSignalNorm /  # noqa: W504
                 curSumm["PSignalNorm"]) if curSumm["PSignalNorm"] != 0
                else Decimal(0))
            curAccession.PSignalAndScNormRatiosAverage = (
                (curAccession.ScNormToFileNormRatio +  # noqa: W504
                 curAccession.PSignalNormToFileNormRatio) / 2)


def GetScPsigAndNormFilesSumm(
    accessionsPerTable: Dict[str, Dict[str, Accession]]
) -> Dict[str, Dict[str, Decimal]]:
    """ Получаем суммы параметров Sc, Sequence, PrecursorSignal,
    ScNorm, PSignalNorm по файлам


    Возвращаем словарь с суммами вида
    {
        "1.1": {
            "ScSumm": S1,
            "PSignalSumm": P1,
            "ScNorm: SN1,
            "PSignalNorm": PN1
        },
        "1.2": {
            "ScSumm": S2,
            "PSignalSumm": P2,
            "ScNorm: SN2,
            "PSignalNorm": PN2
        },
        ...,
        "n-ый файл": {
            "ScSumm": Sn,
            "PSignalSumm": Pn,
            "ScNorm: SNn,
            "PSignalNorm": PNn
        }
    }
    """

    fileSumms: Dict[str, Dict[str, Decimal]] = {}

    for tableNum, curTable in accessionsPerTable.items():
        fileSumms[tableNum] = {}
        curSumm = fileSumms[tableNum]
        curSumm["ScSumm"] = Decimal(0)
        curSumm["PSignalSumm"] = Decimal(0)
        curSumm["ScNorm"] = Decimal(0)
        curSumm["PSignalNorm"] = Decimal(0)

        for accession in curTable:
            curSumm["ScSumm"] += Decimal(curTable[accession].ScSumm)
            curSumm["ScNorm"] += Decimal(curTable[accession].ScNorm)
            curSumm["PSignalSumm"] += Decimal(curTable[accession].PSignalSumm)
            curSumm["PSignalNorm"] += Decimal(curTable[accession].PSignalNorm)
    return fileSumms


def ApplyWhiteList(peptideTables: Dict[str, Dict[str, List[str]]],
                   whiteList: List[str],
                   columnNames: ColumnNames) -> None:
    """ Удаляем из всех таблиц все id, отсутствующие в белом списке """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable[columnNames.accession])
        i = 0
        while i < curTableLen:
            if (curTable[columnNames.accession][i].split(';')[0] not in
                    whiteList):
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1


def RemoveAccessionsFromTableByBlacklist(peptideTable: Dict[str, List[str]],
                                         blackList: List[str],
                                         columnNames: ColumnNames) -> None:
    curTableLen = len(peptideTable[columnNames.accession])
    i = 0
    while i < curTableLen:
        if(peptideTable[columnNames.accession][i].split(';')[0] in blackList):
            RemoveRow(peptideTable, i)
            curTableLen -= 1
            i -= 1
        i += 1


def ApplyBlackList(peptideTables: Dict[str, Dict[str, List[str]]],
                   blackList: List[str],
                   columnNames: ColumnNames) -> None:
    """ Удаляем из всех таблиц все id, находящиеся в чёрном списке """

    for curTable in peptideTables.values():
        RemoveAccessionsFromTableByBlacklist(
            curTable, blackList, columnNames)


def TestUnusedContribConfParams(unused: Comparable,
                                contrib: Comparable,
                                conf: Comparable,
                                peptideTable: Dict[str, List[str]],
                                tableNum: str,
                                tableRowNum: int,
                                columnNames: ColumnNames) -> bool:
    if(unused.compare(
        peptideTable[columnNames.unused][tableRowNum], tableNum) and  # noqa: W504, E501
       contrib.compare(
           peptideTable[columnNames.contribution][tableRowNum], tableNum) and  # noqa: W504, E501
       conf.compare(
           peptideTable[columnNames.confidence][tableRowNum], tableNum)):
        return True
    return False


def ApplyParamsFilter(unused: Comparable,
                      contrib: Comparable,
                      conf: Comparable,
                      peptideTables: Dict[str, Dict[str, List[str]]],
                      columnNames: ColumnNames) -> None:
    """ Применяем фильтры, завязанные на параметры unused, contib, conf.
    """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable[columnNames.unused])
        i = 0

        while i < curTableLen:
            if not TestUnusedContribConfParams(
                    unused, contrib, conf, curTable, tableNum, i, columnNames):
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1


def TestConfDefaultCondition(confVal: Decimal) -> int:
    if confVal >= Decimal("99"):
        return 2
    if confVal >= Decimal("95"):
        return 1
    return 0


## Создание словаря, содержащего все accession из данной таблицы в качестве
# ключей и 0 в качестве значений
# @param peptideTable Таблица-результат чтения PeptideSummary файла с помощью
#                     Classes.ReadTable
# @param columnNames Экземпляр класса Classes.ColumnNames
# @returns Словарь, содержащий все accession данной таблицы в качестве ключей
#          и 0 в качестве значений
def GenerateTableAccessionsBunch(peptideTable: Dict[str, List[str]],
                                 columnNames: ColumnNames) -> Dict[str, int]:
    accessionBunch: Dict[str, int] = {}
    for accession in peptideTable[columnNames.accession]:
        if accession not in accessionBunch:
            accessionBunch[accession] = 0
    return accessionBunch


def ApplyConfidenceDefaultFilter(
        peptideTables: Dict[str, Dict[str, List[str]]],
        columnNames: ColumnNames) -> None:
    """ Применяем Confidence default условие.

    Удаляем все ID, у которых нет ни одной строки проходящей default-условие.

    default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95"""
    for tableNum, curTable in peptideTables.items():
        curTableLen = len(curTable[columnNames.unused])
        # Заносим все Accession в список Accession для удаления
        # В процессе чтения списка записи из него будут удаляться
        blackList = GenerateTableAccessionsBunch(curTable, columnNames)
        i = 0
        while i < curTableLen:
            curAccession = curTable[columnNames.accession][i]
            if curAccession in blackList:
                blackList[curAccession] += (
                    TestConfDefaultCondition(
                        Decimal(curTable[columnNames.confidence][i])))
                if blackList[curAccession] > 1:
                    del blackList[curAccession]
            i += 1
        RemoveAccessionsFromTableByBlacklist(curTable,
                                             [*blackList],
                                             columnNames)


def ApplyConfidenceIDFilter(confID: Comparable,
                            peptideTables: Dict[str, Dict[str, List[str]]],
                            columnNames: ColumnNames) -> None:
    """ Применяем Confidence ID фильтр.

    Удаляем все ID, у которых нет ни одной строки проходящей условия,
    заданные confID.

    default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95"""
    if confID.val is None:
        ApplyConfidenceDefaultFilter(peptideTables, columnNames)
    else:
        for tableNum, curTable in peptideTables.items():
            curTableLen = len(curTable[columnNames.unused])
            # Заносим все Accession в список Accession для удаления
            # В процессе чтения списка записи из него будут удаляться
            blackList = GenerateTableAccessionsBunch(curTable, columnNames)
            i = 0
            while i < curTableLen:
                curAccession = curTable[columnNames.accession][i]
                if curAccession in blackList:
                    if confID.compare(curTable[columnNames.confidence][i],
                                      tableNum):
                        del blackList[curAccession]
                i += 1
            RemoveAccessionsFromTableByBlacklist(curTable,
                                                 [*blackList],
                                                 columnNames)


def GetFileLines(filename: str) -> Union[List[str], None]:
    """ Returns list of file strings without newline symbols """

    if len(filename):
        with open(filename) as tfile:
            return tfile.read().split('\n')
    return None


def ReadSeqDB(seqDBFilename: str) -> Dict[str, Sequence]:
    """ Считывание последовательностей из файла

    Считывание длин последовательностей из файла БД с последовательностями в
    словарь классов Sequence вида {"Accession": Sequence} """

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
                    while ((i < len(strings)) and (
                            not len(strings[i]) or strings[i][0] != '>')):
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


def GetInput() -> Input:
    inputParams = Input()
    inputParams.inputPath = "Input"
    if len(argv) == 12:
        inputParams.proteinPilotVersion = argv[1]
        inputParams.whiteList = GetFileLines(argv[2])
        inputParams.blackList = GetFileLines(argv[3])
        inputParams.isProteinGroupFilter = argv[4].strip().lower()
        inputParams.seqDB = ReadSeqDB(argv[5])
        inputParams.unused = Comparable(argv[6])
        inputParams.contrib = Comparable(argv[7])
        inputParams.confID = argv[8]
        inputParams.confPeptide = argv[9]
        inputParams.minGroupsWithAccession = int(argv[10])
        inputParams.maxGroupAbsence = int(argv[11])
    else:
        inputParams.proteinPilotVersion = input(
            "ProteinPilot Version (4 or 5): ")
        inputParams.whiteList = GetFileLines(input("ID list file name: "))
        inputParams.blackList = GetFileLines(
            input("ID exclusion list file name: "))
        inputParams.isProteinGroupFilter = input(
            "Protein group filter (Y or N): ").strip()
        inputParams.seqDB = ReadSeqDB(input("Database file name: "))
        inputParams.unused = Comparable(input("Unused: "))
        inputParams.contrib = Comparable(input("Contribution: "))
        inputParams.confID = input("Confidence ID: ")
        inputParams.confPeptide = input("Confidence peptide: ")
        inputParams.minGroupsWithAccession = int(input("Min groups with ID: "))
        inputParams.maxGroupAbsence = int(
            input("Max missing values per group: "))
    return inputParams
