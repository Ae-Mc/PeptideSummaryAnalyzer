from dataclasses import dataclass
from decimal import Decimal


@dataclass
class Accession:
    """Хранит информацию об Accession

    Attributes:
        name: имя Accession
        Unused: параметр Unused Accession
        ScSumm: сумма параметра Sc для всех строк с данным Accession в файле
        ScNorm: отношение ScSumm к длинне последовательности данного Accession
            в базе данных последовательностей
        ScNormToFileNormRatio: отношение ScNorm к сумме всех ScNorm в файле
        PSignalSumm: сумма параметра Precursor Signal для всех строк с данным
            Accession в файле
        PSignalNorm: отношение PSignalSumm к длинне последовательности данного
            Accession в базе данных последовательностей
        PSignalNormToFileNormRatio: отношение PSignalNorm к сумме всех
            PSignalNorm в файле
        SeqlenSumm: сумма параметра Sequence для всех строк с данным
            Accession в файле
        Counts: количество строк с данным Accession в файле
        PSignalAndScNormRatiosAverage: среднее арифметическое
            ScNormToFileNormRatio и PSignalNormToFileNormRatio
    """
    name: str = ""
    ScSumm: Decimal = Decimal(0)
    ScNorm: Decimal = Decimal(0)
    ScNormToFileNormRatio: Decimal = Decimal(0)
    PSignalSumm: Decimal = Decimal(0)
    PSignalNorm: Decimal = Decimal(0)
    PSignalNormToFileNormRatio: Decimal = Decimal(0)
    SeqlenSumm: int = 0
    Counts: int = 0
    PSignalAndScNormRatiosAverage: Decimal = Decimal(0)
