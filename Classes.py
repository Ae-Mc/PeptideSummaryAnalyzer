from typing import Dict, List


def IsFloat(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


# Класс для хранения последовательностей
class Sequence:
    len: int
    desc: str
    __seq: str

    @property
    def seq(self):
        return self.__seq

    @seq.setter
    def seq(self, seq: str):
        self.len = len(seq)
        self.__seq = seq

    def __init__(self) -> None:
        self.len = 0
        self.desc = ""
        self.seq = ""


# Класс для хранения параметров фильтра
class Comparable:
    op: str
    val: str
    __float: bool

    def __init__(self, op: str = "", unused: str = "") -> None:
        self.op: str = op
        self.val: str = unused
        self.__float: bool = False

    @property
    def val(self):
        return self.__val

    @val.setter
    def val(self, value):
        self.__val = value
        try:
            float(value)
            self.__float = True
        except (TypeError, ValueError):
            self.__float = False

    def GetComparable(self, paramString: str):
        """ Получение param

        Получение param, общего для всех фалов, либо для каждого своего, из
        строки.
        Формат: [операция][[имя файла] или [число]]
        Примеры:
            >=99
            < paramList.txt"""

        paramString = paramString.strip()
        self.op = ''.join([ch for ch in paramString if ch in "!=<>"])

        paramString = paramString[len(self.op):].strip()
        if IsFloat(paramString):
            self.val = float(paramString)
        elif len(self.op) and len(paramString):
            with open(paramString) as paramStringFile:
                strings = paramStringFile.read().replace(' ', '\t').split('\n')
                paramStrings = {}
                for string in strings:
                    string = string.strip()
                    if len(string):
                        paramStrings[string.split('\t')[0]] = (
                            string.split('\t')[1])
                self.val = paramStrings

    @staticmethod
    def GetComparableClass(paramString: str):
        """ Получение param

        Получение param, общего для всех фалов, либо для каждого своего, из
        строки.
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
                        paramStrings[string.split('\t')[0]] = (
                            string.split('\t')[1])
                param.val = paramStrings
        return param

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


# Класс для хранения параметров Accession
class Accession:
    Unused: float
    ScSumm: float
    ScNorm: float
    ScNormToFileNormRatio: float
    PSignalSumm: float
    PSignalNorm: float
    PSignalNormToFileNormRatio: float
    SeqlenSumm: int
    Counts: int
    PSignalAndScNormRatiosAverage: float

    def __init__(self, name="", Unused=0, ScSumm=0, ScNorm=0,
                 ScNormToFileNormRatio=0, PSignalSumm=0, PSignalNorm=0,
                 PSignalNormToFileNormRatio=0, SeqlenSumm=0, Counts=0,
                 PSignalAndScNormRatiosAverage=0):
        self.name = name
        self.Unused = Unused
        self.ScSumm = ScSumm
        self.ScNorm = ScNorm
        self.ScNormToFileNormRatio = ScNormToFileNormRatio
        self.PSignalSumm = PSignalSumm
        self.PSignalNorm = PSignalNorm
        self.PSignalNormToFileNormRatio = PSignalNormToFileNormRatio
        self.SeqlenSumm = SeqlenSumm
        self.Counts = Counts
        self.PSignalAndScNormRatiosAverage = PSignalAndScNormRatiosAverage


class Input:
    seqDB: Dict[str, Sequence]
    unused: Comparable
    contrib: Comparable
    __conf: Comparable
    isDefaultConf: bool
    whiteList: List[str]
    blackList: List[str]
    minGroupsWithAccession: int
    maxGroupAbsence: int

    @property
    def conf(self):
        return self.__conf

    @conf.setter
    def conf(self, val):
        self.__conf = Comparable.GetComparableClass(val.strip())
        if val.strip() == "":
            self.isDefaultConf = True
        else:
            self.isDefaultConf = False
