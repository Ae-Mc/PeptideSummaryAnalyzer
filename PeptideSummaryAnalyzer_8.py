from typing import List, Dict, Union
from os import listdir

INPUTPATH = "./Input"


# Класс для хранения последовательностей
class Sequence():
    len: int
    desc: str
    seq: str

    def __init__(self) -> None:
        self.len = 0
        self.desc = ""
        self.seq = ""


# Класс для хранения параметров фильтра
class Comparable():
    op: str
    val: str
    __float: bool

    def __init__(self, op: str="", unused: str="") -> None:
        self.op: str = op
        self.val: str = unused
        self.__float: bool = False

    @property
    def val(self):
        return self.__val

    @val.setter
    def val(self, value):
        self.__val = value
        self.__float = IsFloat(value)

    def compare(self, value, filename: str) -> bool:
        if self.__float:
            return eval("{value}{op}{comparable}".format(
                value=value,
                op=self.op,
                comparable=self.val))
        else:
            if len(self.val) and filename in self.val:
                return eval("{value}{op}{comparable}".format(
                    value=value,
                    op=self.op,
                    comparable=self.val[filename]))
            else:
                return True


def ReadSeqDB(seqDBFilename: str) -> dict:
    """ Считывание последовательностей из файла

    Считывание длин последовательностей из файла БД с последовательностями в
    словарь массивов вида {"id"} = {"seqLen", "seq", "seqDesc}"""

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
            i += 1

        for seqID in seqDB:
            if not seqDB[seqID].len:
                input("Error!\nLength of sequence with id " + seqID + " = 0")
                raise(IndexError)
        return seqDB
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


def GetFileLines(filename: str) -> list:
    """ Returns list of file strings without newline symbols """

    try:
        with open(filename) as tfile:
            return tfile.read().split('\n')
    except FileNotFoundError:
        return None


def ApplyDefaultConf(peptideTables: Dict[str, Dict[str, List[str]]]):
    """ default-условие: Accession учитывается, если хотя бы одна
    из строк с ним имеет Conf >= 99 или минимум две строки имеют
    Conf >= 95"""

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Unused"])
        # Заносим все Accession в список Accession для удаления
        # В процессе чтения списка записи из него будут удаляться
        blackListByConfCondition: Dict[str, int] = {}
        for accession in curTable["Accessions"]:
            if accession not in blackListByConfCondition:
                blackListByConfCondition[
                    accession.split(';')[0]] = 0

        i = 0
        while i < curTableLen:
            curAccession = curTable["Accessions"][i].split(';')[0]
            if curAccession in blackListByConfCondition:
                    blackListByConfCondition[curAccession] += (
                        TestConfDedaultCondition(curTable["Conf"][i]))
                    if blackListByConfCondition[curAccession] > 1:
                        del blackListByConfCondition[curAccession]
        RemoveAccessionsFromTableByBlacklist(curTable, blackListByCondition)
    return peptideTables

def RemoveAccessionsFromTableByBlackList(peptideTable: Dict[str, List[str]],
    blackList: Dict[str, List[str]]):
    i = 0
    while i < curTableLen:
        if(peptideTable["Accessions"][i].split(';')[0] in blackList):
            RemoveRow(curTable, i)
            curTableLen -= 1
            i -= 1
        i += 1


def TestConfDefaultCondition(confVal: int):
    if float(confVal) >= 99.0:
        return 2
    if float(confVal) >= 95.0:
        return 1
    return 0


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


def ApplyParamsFilter(unused: Comparable,
                      contrib: Comparable,
                      conf: Comparable,
                      peptideTables: Dict[str, Dict[str, List[str]]]) -> dict:
    """ Применяем фильтры, завязанные на параметры unused, contib, conf.

    peptideTables - словарь вида
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

    Возвращаем словарь с уже применёнными фильрами.
    """

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        curTableLen = len(curTable["Unused"])
        i = 0

        while i < curTableLen:
            if not (unused.compare(curTable["Unused"][i], tableNum) and
                    contrib.compare(curTable["Contrib"][i], tableNum) and
                    conf.compare(curTable["Conf"][i], tableNum)):
                RemoveRow(curTable, i)
                i -= 1
                curTableLen -= 1
            i += 1

    return peptideTables


# TODO: Сделать применение белого/чёрного списков
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


# TODO: сделать подсчёт Sc, Seq и Psignal сумм
def GetScPsigFileSumm(peptideTables: Dict[str, Dict[str, List[str]]]) -> Dict:
    """ Получаем суммы параметров Sc, Sequence, PrecursorSignal по файлам

    peptideTables - словарь вида
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

    Возвращаем словарь с суммами вида
    {
        "1.1": {
            "Sc": S1,
            "pSignal": P1
        },
        "1.2": {
            "Sc": S2,
            "pSignal": P2
        },
        ...,
        "N-ый файл": {
            "Sc": SN,
            "pSignal": PN
        }
    }
    """

    fileSumms: Dict[str, Dict[str, Union[float, int]]] = {}

    for tableNum in peptideTables:
        curTable = peptideTables[tableNum]
        fileSumms[tableNum] = {}
        curSumm = fileSumms[tableNum]
        curSumm["Sc"] = 0
        curSumm["pSignal"] = 0

        i = 0
        while i < len(curTable["Sc"]):
            curSumm["Sc"] += float(curTable["Sc"][i])
            curSumm["pSignal"] += float(curTable["PrecursorSignal"][i])
    return fileSumms

# Может ли переменная быть приведена к типу float
def IsFloat(var: str) -> bool:
    try:
        float(var)
        return True
    except (ValueError, TypeError):
        return False


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
    unused = GetParam(">=unusedList.txt")
    contrib = GetParam(">=1")
    conf = ">=50"
    whiteList = GetFileLines("")
    blackList = GetFileLines("IDexcl.txt")

    if conf.strip() == "default":
        ApplyDefaultConf(peptideTables)
        ApplyParamsFilter(unused, contrib, GetParam(""), peptideTables)
    else:
        ApplyParamsFilter(unused, contrib, GetParam(conf), peptideTables)

    # minGroupsWithId = 20
    # maxGroupAbsence = 5

    ApplyBlackList(peptideTables, blackList)
    fileSumms = GetScPsigFileSumm(peptideTables)


if __name__ == "__main__":
    main()
