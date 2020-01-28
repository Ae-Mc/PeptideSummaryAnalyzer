# Класс для хранения последовательностей
class Sequence:
    len: int
    desc: str
    seq: str

    def __init__(self) -> None:
        self.len = 0
        self.desc = ""
        self.seq = ""


# Класс для хранения параметров фильтра
class Comparable:
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
        try:
            float(value)
            self.__float = True
        except (TypeError, ValueError):
            self.__float = False

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
    unused: float
    ScSumm: float
    ScNorm: float
    ScNormToFileNormRatio: float
    PSignalSumm: float
    PSignalNorm: float
    PSignalNormToFileNormRatio: float
    SeqlenSumm: int
    Counts: int
    PSignalAndScNormRatiosAverage: float

    def __init__(self, name="", unused=0, ScSumm=0, ScNorm=0,
                 ScNormToFileNormRatio=0, PSignalSumm=0, PSignalNorm=0,
                 PSignalNormToFileNormRatio=0, SeqlenSumm=0, Counts=0,
                 PSignalAndScNormRatiosAverage=0):
        self.name = name
        self.unused = unused
        self.ScSumm = ScSumm
        self.ScNorm = ScNorm
        self.ScNormToFileNormRatio = ScNormToFileNormRatio
        self.PSignalSumm = PSignalSumm
        self.PSignalNorm = PSignalNorm
        self.PSignalNormToFileNormRatio = PSignalNormToFileNormRatio
        self.SeqlenSumm = SeqlenSumm
        self.Counts = Counts
        self.PSignalAndScNormRatiosAverage = PSignalAndScNormRatiosAverage
