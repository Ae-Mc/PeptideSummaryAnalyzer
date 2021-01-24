from sys import argv
from typing import List, Dict, Union
from .Accession import Accession
from .AccessionTables import AccessionTables
from .Comparable import Comparable
from .Input import Input
from .Sequence import Sequence
from .PeptideTable import PeptideTable
from .PeptideTables import PeptideTables
from decimal import Decimal


def RemoveRow(table: Dict[str, List[str]], rowNum: int) -> None:
    """Удаляет строку из таблицы, полученной с помощью .ReadTable.ReadTable

    Args:
        table Таблица, из которой происходит удаление
        rowNum Номер строки в таблице
    """
    columns = [column for column in table]
    for column in columns:
        del table[column][rowNum]


def CountAccessionLackInGroup(
        accession: str,
        group: List[str],
        accessionsPerTable: Dict[str, Dict[str, Accession]]) -> int:
    """Подсчитывает количество файлов, не содержащих данный accession, в группе

    Args:
        accession: Имя accession, для которого считается отсутствие в группе
        group: Список номеров таблиц, входящих в группу
        accessionsPerTable: Словарь, ключом в котором является номер таблицы, а
            ключом — словарь, в котором, в свою очередь, ключом
            является имя accession, а значением — экземпляр
            класса Classes.Accession

    Returns:
        Количество файлов, не содержащих данный accession, в группе
    """
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
    """Подсчёт количества групп, в которых данный accession отсутствует не
    более maxGroupAbsence раз

    Args:
        groups: Словарь, в котором ключом является номер группы, а значением —
                  список номеров входящих в неё таблиц
        accession: Имя accession, для которого считается отсутствие в группах
        maxGroupAbsence: Максимальное количество отсутствий accession в группе
        accessionsPerTable: Словарь, ключом в котором является номер таблицы, а
            ключом — словарь, в котором, в свою очередь, ключом является имя
            accession, а значением — экземпляр Accession

    Returns:
        Количество групп, в которых accession отсутствует не более
        maxGroupAbsence раз
    """

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
    """Создаёт словарь с группами, в котором ключом является номер группы, а
    значением — список номеров таблиц, входящих в группу

    Args:
        accessionsPerTable: Словарь, ключом в котором является номер таблицы, а
            ключом — словарь, в котором, в свою очередь, ключом
            является имя accession, а значением — экземпляр
            класса Classes.Accession

    Returns:
        Словарь с группами, в котором ключом является номер группы, а
        значением — список номеров таблиц, входящих в группу
    """

    for tableNum in accessionsPerTable.keys():
        groupNum = tableNum.split('.')[0]
        if groupNum not in groups:
            groups[groupNum] = []
        groups[groupNum].append(tableNum)
    return groups


def ApplyGroupFilter(accessionTables: AccessionTables,
                     maxGroupAbsence: int,
                     minGroupsWithAccession: int) -> None:
    """Применение фильтра по группам к AccessionTables

    Таблицы, начинающиеся с одинакового числа входят в одну группу, например:
    1.1, 1.2, 1.3 входят в одну группу, 2.1, 2.2 входят в другую и т. д.
    Из всех таблиц удаляются все Accession, присутствующие в меньше, чем в
    minGroupsWithAccession группах. Accession считается отсутствуюющим в
    группе, если он отсутствует в больше, чем maxGroupAbsence таблицах внутри
    группы.

    Args:
        accessionTables: Экземпляр класса Classes.AccessionTables, к которому
            применяется фильтр
        maxGroupAbsence: максимально возможное количество таблиц без Accession
            в группе
        minGroupsWithAccession: минимальное количество групп с accession
    """
    groups: Dict[str, List[str]] = GenerateGroupsBunch(accessionTables)

    accessionBunch = accessionTables.GenerateAccessionsBunchOverAllTables()
    for accession in accessionBunch:
        if CountGroupsWithAccession(
                groups,
                accession,
                maxGroupAbsence,
                accessionTables) < minGroupsWithAccession:
            accessionTables.RemoveAccessionFromAllTables(accession)


def CalculateAccessionsNormRatios(
        accessionTables: AccessionTables,
        tableSumms: Dict[str, Dict[str, Decimal]]) -> None:
    """Подсчёт нормализованных значений PSignalNormToFileNormRatio и

    Нормализованные значения PSignalNormToFileNormRatio и ScNormToFileNormRatio
    это отношения ScNorm/PsignalNorm accession к сумме всех ScNorm/PsignalNorm
    в файле ScNormToFileNormRatio для всех accession

    Args:
        accessionTables: Экземпляр класса Classes.AccessionTables, к которому
            применяется фильтр
        tableSumms: Словарь, содержащий суммы значений ScNorm и PSignalNorm
            для каждой таблицы
    """
    for tableNum, curSumm in tableSumms.items():
        for accession in accessionTables[tableNum].values():
            accession.ScNormToFileNormRatio = (
                (accession.ScNorm
                 / curSumm["ScNorm"]) if curSumm["ScNorm"] != 0
                else Decimal(0))
            accession.PSignalNormToFileNormRatio = (
                (accession.PSignalNorm
                 / curSumm["PSignalNorm"]) if curSumm["PSignalNorm"] != 0
                else Decimal(0))
            accession.PSignalAndScNormRatiosAverage = (
                (accession.ScNormToFileNormRatio
                 + accession.PSignalNormToFileNormRatio) / 2)


def GetScPsigAndNormFilesSumm(
    accessionsPerTable: Dict[str, Dict[str, Accession]]
) -> Dict[str, Dict[str, Decimal]]:
    """Получает суммы параметров Sc, Sequence, PrecursorSignal,
    ScNorm, PSignalNorm по файлам

    Args:
        accessionsPerTable: Словарь {
                "Номер таблицы": {
                    "Имя Accession": Accession
                }
            }

    Returns:
        Словарь с суммами {
            "Номер таблицы": {
                "ScSumm": x1,
                "PSignalSumm": x2,
                "ScNorm: x3,
                "PSignalNorm": x4
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


def RemoveAccessionsListFromTable(peptideTable: PeptideTable,
                                  blackList: List[str]) -> None:
    """Удаляет из таблицы все id, находящиеся в списке blackList

    Args:
        peptideTable: Peptide таблица
        blackList: список Accession для удаления
    """
    i = 0
    while i < len(peptideTable):
        if(peptideTable[i].name.split(';')[0] in blackList):
            peptideTable.pop(i)
            continue
        i += 1


def ApplyBlackList(peptideTables: PeptideTables,
                   blackList: List[str]) -> None:
    """Удаляет из всех таблиц все id, находящиеся в чёрном списке

    Args:
        peptideTables: словарь с таблицами Peptide вида: {
                "номер таблицы": PeptideTable
            }
        blackList: список Accession, находящихся в чёрном списке
    """

    for curTable in peptideTables.values():
        RemoveAccessionsListFromTable(curTable, blackList)


def TestUnusedContribConfParams(unused: Comparable,
                                contrib: Comparable,
                                conf: Comparable,
                                peptideTable: PeptideTable,
                                tableNum: str,
                                tableRowNum: int) -> bool:
    """Проверка unused, contribution и confidence параметров

    Args:
        unused: параметр фильтра unused
        contrib: параметр фильтра contribution
        conf: параметр фильтра confidence
        peptideTable: Peptide таблица
        tableNum: номер таблицы, в которой находится значение
        tableRowNum: номер строк в таблице
    """
    if(unused.compare(
        peptideTable[tableRowNum].unused, tableNum)
       and contrib.compare(
           peptideTable[tableRowNum].contribution, tableNum)
       and conf.compare(
           peptideTable[tableRowNum].confidence, tableNum)):
        return True
    return False


def ApplyParamsFilter(unused: Comparable,
                      contrib: Comparable,
                      conf: Comparable,
                      peptideTables: PeptideTables) -> None:
    """Применяем фильтры, завязанные на параметры unused, contribution,
    confidence

    Args:
        unused: параметр фильтра unused
        contrib: параметр фильтра contribution
        conf: параметр фильтра confidence
        peptideTables: словарь с таблицами Peptide вида: {
                "номер таблицы": PeptideTable
            }
    """

    for tableNum, curTable in peptideTables.items():
        i = 0

        while i < len(curTable):
            if not TestUnusedContribConfParams(
                    unused, contrib, conf, curTable, tableNum, i):
                curTable.pop(i)
                continue
            i += 1


def TestConfDefaultCondition(confVal: Decimal) -> int:
    """Проверка default-условия для confidence

    Args:
        confVal: значение confidence, которое проверяется

    Returns:
        2, если confVal >= 99
        1, если confVal >= 95
        иначе 0
    """
    if confVal >= Decimal("99"):
        return 2
    if confVal >= Decimal("95"):
        return 1
    return 0


def GenerateTableAccessionsBunch(peptideTable: PeptideTable) -> Dict[str, int]:
    """Создание словаря, содержащего все accession из данной таблицы в качестве
    ключей и 0 в качестве значений

    Args:
        peptideTable: Таблица-результат чтения PeptideSummary файла с помощью
            Classes.ReadTable

    Returns:
        Словарь, содержащий все accession данной таблицы в качестве ключей
        и 0 в качестве значений
    """
    accessionBunch: Dict[str, int] = {}
    for line in peptideTable:
        if line.name not in accessionBunch:
            accessionBunch[line.name] = 0
    return accessionBunch


def ApplyConfidenceDefaultFilter(peptideTables: PeptideTables) -> None:
    """ Применяем Confidence default условие.

    Удаляем все ID, у которых нет ни одной строки проходящей default-условие.

    default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95

    Args:
        peptideTables: словарь вида: {
                "Номер таблицы": PeptideTable
            }
    """
    for tableNum, curTable in peptideTables.items():
        curTableLen = len(curTable)
        # Заносим все Accession в список Accession для удаления
        # В процессе чтения списка записи из него будут удаляться
        blackList = GenerateTableAccessionsBunch(curTable)
        i = 0
        while i < curTableLen:
            curAccession = curTable[i].name
            if curAccession in blackList:
                blackList[curAccession] += (
                    TestConfDefaultCondition(
                        Decimal(curTable[i].confidence)))
                if blackList[curAccession] > 1:
                    del blackList[curAccession]
            i += 1
        RemoveAccessionsListFromTable(curTable, [*blackList])


def ApplyConfidenceIDFilter(confID: Comparable,
                            peptideTables: PeptideTables) -> None:
    """ Применяем Confidence ID фильтр.

    Удаляем все ID, у которых нет ни одной строки проходящей условия,
    заданные confID.

    default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95

    Args:
        confID: параметр фильтра confidence
        peptideTables: словарь с таблицами Peptide вида: {
                "номер таблицы": PeptideTables
            }
    """
    if confID.val is None:
        ApplyConfidenceDefaultFilter(peptideTables)
    else:
        for tableNum, curTable in peptideTables.items():
            curTableLen = len(curTable)
            # Заносим все Accession в список Accession для удаления
            # В процессе чтения списка записи из него будут удаляться
            blackList = GenerateTableAccessionsBunch(curTable)
            i = 0
            while i < curTableLen:
                curAccession = curTable[i].name
                if curAccession in blackList:
                    if confID.compare(curTable[i].confidence, tableNum):
                        del blackList[curAccession]
                i += 1
            RemoveAccessionsListFromTable(curTable, [*blackList])


def GetFileLines(filename: str) -> Union[List[str], None]:
    """Считывает файл и возвращает список строк без символа переноса строки

    Args:
        filename: Имя файла, из которого происходит считывание строк

    Returns:
        Список строк файла, без символов переноса строки
    """

    if len(filename):
        with open(filename) as tfile:
            return tfile.read().split('\n')
    return None


def ReadSeqDB(seqDBFilename: str) -> Dict[str, Sequence]:
    """ Считывание последовательностей из файла

    Считывание длин последовательностей из файла БД с последовательностями в
    словарь классов Sequence вида {"Accession": Sequence}

    Args:
        seqDBFilename: Имя файла с последовательностями
    Returns:
        Словарь вида {"Accession": Sequence}
    """

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
    """Получение параметров для запуска обработки

    Returns:
        Класс Input, содержащий все нужные параметры для запуска обработки
    """
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
