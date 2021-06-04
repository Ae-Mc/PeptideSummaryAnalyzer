from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTable import RawPeptideTable
from Classes.RawPeptideTables import RawPeptideTables
from decimal import Decimal
from os import listdir, path
from sys import argv
from typing import List, Dict, Union, Optional, Set
from .Accession import Accession
from .AccessionTables import AccessionTables
from .Comparable import Comparable
from .Input import Input
from .ProteinAccession import ProteinAccession
from .ProteinTable import ProteinTable
from .PeptideAccession import PeptideAccession
from .PeptideTable import PeptideTable
from .PeptideTables import PeptideTables
from .SequenceDatabase import SequenceDatabase


def CountAccessionLackInGroup(
    accession: str,
    group: List[str],
    accessionsPerTable: Dict[str, Dict[str, Accession]],
) -> int:
    """Подсчитывает количество файлов, не содержащих данный accession, в группе

    Args:
        accession: Имя accession, для которого считается отсутствие в группе
        group: Список номеров таблиц, входящих в группу
        accessionsPerTable: Словарь, ключом в котором является номер таблицы, а
            значением — словарь, в котором, в свою очередь, ключом
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
    accessionsPerTable: Dict[str, Dict[str, Accession]],
) -> int:
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
        if (
            CountAccessionLackInGroup(
                accession, tableNames, accessionsPerTable
            )
            > maxGroupAbsence
        ):
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
        groupNum = tableNum.split(".")[0]
        if groupNum not in groups:
            groups[groupNum] = []
        groups[groupNum].append(tableNum)
    return groups


def ApplyGroupFilter(
    accessionTables: AccessionTables,
    maxGroupAbsence: int,
    minGroupsWithAccession: int,
) -> None:
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
        if (
            CountGroupsWithAccession(
                groups, accession, maxGroupAbsence, accessionTables
            )
            < minGroupsWithAccession
        ):
            accessionTables.RemoveAccessionFromAllTables(accession)


def CalculateAccessionsNormRatios(
    accessionTables: AccessionTables, tableSumms: Dict[str, Dict[str, Decimal]]
) -> None:
    """Подсчёт нормализованных значений PSignalNormToFileNormRatio и
    ScNormToFileNormRatio

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
                (accession.ScNorm / curSumm["ScNorm"])
                if curSumm["ScNorm"] != 0
                else Decimal(0)
            )
            accession.PSignalNormToFileNormRatio = (
                (accession.PSignalNorm / curSumm["PSignalNorm"])
                if curSumm["PSignalNorm"] != 0
                else Decimal(0)
            )


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


def RemoveAccessionsListFromTable(
    peptideTable: PeptideTable, blackList: List[str]
) -> None:
    """Удаляет из таблицы все id, находящиеся в списке blackList

    Args:
        peptideTable: Peptide таблица
        blackList: список Accession для удаления
    """
    i = 0
    while i < len(peptideTable):
        if peptideTable[i].name.split(";")[0] in blackList:
            peptideTable.pop(i)
            continue
        i += 1


def ApplyBlackList(
    rawPeptideTables: RawPeptideTables, blackList: List[str]
) -> None:
    """Удаляет из всех таблиц все id, находящиеся в чёрном списке

    Args:
        rawPeptideTables: словарь с таблицами Peptide вида: {
                "номер таблицы": RawPeptideTable
            }
        blackList: список Accession, находящихся в чёрном списке
    """

    curTable: RawPeptideTable
    for curTable in rawPeptideTables.values():
        curTable.RemoveRowsWithAccessions(blackList)


def TestConfParams(
    conf: Comparable,
    peptideTable: PeptideTable,
    tableNum: str,
    tableRowNum: int,
) -> bool:
    """Проверка confidence параметра

    Args:
        conf: параметр фильтра confidence
        peptideTable: Peptide таблица
        tableNum: номер таблицы, в которой находится значение
        tableRowNum: номер строк в таблице
    """
    if conf.compare(peptideTable[tableRowNum].confidence, tableNum):
        return True
    return False


def ApplyPeptideConfidenceFilter(
    conf: Comparable, peptideTables: PeptideTables
) -> None:
    """Применяем фильтр, завязанный на параметр confidence

    Args:
        conf: параметр фильтра confidence
        peptideTables: словарь с таблицами Peptide вида: {
                "номер таблицы": PeptideTable
            }
    """

    for tableNum, curTable in peptideTables.items():
        i = 0

        while i < len(curTable):
            if not TestConfParams(conf, curTable, tableNum, i):
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
    """Применяем Confidence default условие.

    Удаляем все ID, у которых нет ни одной строки проходящей default-условие.

    default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95

    Args:
        peptideTables: словарь вида: {
                "Номер таблицы": PeptideTable
            }
    """
    for curTable in peptideTables.values():
        curTableLen = len(curTable)
        # Заносим все Accession в список Accession для удаления
        # В процессе чтения списка записи из него будут удаляться
        blackList = GenerateTableAccessionsBunch(curTable)
        i = 0
        while i < curTableLen:
            curAccession = curTable[i].name
            if curAccession in blackList:
                blackList[curAccession] += TestConfDefaultCondition(
                    Decimal(curTable[i].confidence)
                )
                if blackList[curAccession] > 1:
                    del blackList[curAccession]
            i += 1
        RemoveAccessionsListFromTable(curTable, [*blackList])


def ApplyProteinConfidenceFilter(
    confID: Comparable, peptideTables: PeptideTables
) -> None:
    """Применяем Confidence ID фильтр.

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

    if not (
        filename.endswith("\\")
        or filename.endswith("/")
        or len(filename.strip()) == 0
    ):
        with open(filename) as tfile:
            return tfile.read().split("\n")
    return None


def FindFastaFile(inputPath: str) -> str:
    """Ищет fasta файл в папке inputPath

    Args:
        inputPath: путь к папке, в которой будет производиться поиск fasta
            файла

    Returns:
        Путь к fasta файлу, составленный из inputPath и имени fasta файла
    """
    for file in listdir(inputPath):
        if file.endswith(".fasta"):
            return path.join(inputPath, file)
    raise FileNotFoundError(f'Fasta file not found in folder "{inputPath}"!')


def TestFastaAccessions(
    seqDB: SequenceDatabase,
    rawPeptideTables: RawPeptideTables,
    proteinTables: Optional[Dict[str, ProteinTable]] = None,
) -> None:
    rawPeptideTable: RawPeptideTable
    notFoundAccessions: Set[str] = set()
    for rawPeptideTable in rawPeptideTables.values():
        peptideRow: PeptideRow
        for peptideRow in rawPeptideTable:
            for accession in peptideRow.accessions:
                if accession not in seqDB and not IsReversed(accession):
                    notFoundAccessions.add(accession)
    if proteinTables is not None:
        proteinTable: ProteinTable
        for proteinTable in proteinTables.values():
            proteinAccession: ProteinAccession
            for proteinAccession in proteinTable:
                if proteinAccession.name not in seqDB and not IsReversed(
                    proteinAccession.name
                ):
                    notFoundAccessions.add(proteinAccession.name)
    if len(notFoundAccessions) > 0:
        print(
            f"ERROR! missing {len(notFoundAccessions)} sequences:\n\t"
            + ("\n\t".join(sorted(notFoundAccessions)))
        )
        raise KeyError("ERROR! Missing sequences (check above)!")


def GetInput() -> Input:
    """Получение параметров для запуска обработки

    Returns:
        Класс Input, содержащий все нужные параметры для запуска обработки
    """
    inputParams = Input()
    inputParams.rootPath = "."
    inputParams.inputPath = "./Input"
    inputParams.seqDB = SequenceDatabase.fromFile(
        FindFastaFile(inputParams.rootPath)
    )
    if len(argv) == 8:
        blackListLines = GetFileLines(argv[1])
        inputParams.fdr = argv[1]
        inputParams.blackList = (
            (argv[2], blackListLines) if blackListLines is not None else None
        )
        inputParams.proteinConfidence = argv[4]
        inputParams.proteinGroupingConfidence = argv[5]
        inputParams.confPeptide = argv[6]
        inputParams.minGroupsWithAccession = int(argv[7])
        inputParams.maxGroupAbsence = int(argv[8])
    else:
        print('"ProteinPilot summary analyzer"')
        print("#Protein filter")
        inputParams.fdr = input(
            "Global FDR critical value (<% k or default): "
        )
        blackListFile = input("ID exclusion list: ")
        blackListLines = GetFileLines(blackListFile)
        inputParams.blackList = (
            (blackListFile, blackListLines)
            if blackListLines is not None
            else None
        )
        inputParams.proteinConfidence = input(
            "Peptide confidence (value or default): "
        )
        inputParams.proteinGroupingConfidence = input(
            "Protein grouping (conf): "
        )
        print("#Peptide filter")
        inputParams.confPeptide = input("Peptide confidence (value): ")
        print("#Output filter")
        inputParams.minGroupsWithAccession = int(input("Min groups with ID: "))
        inputParams.maxGroupAbsence = int(
            input("Max missing values per group: ")
        )
    return inputParams


def IsReversed(accession: str) -> bool:
    return accession.startswith("RRRRR")
