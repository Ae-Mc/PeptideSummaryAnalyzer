from decimal import Decimal


class Accession:
    Unused: Decimal
    ScSumm: Decimal
    ScNorm: Decimal
    ScNormToFileNormRatio: Decimal
    PSignalSumm: Decimal
    PSignalNorm: Decimal
    PSignalNormToFileNormRatio: Decimal
    SeqlenSumm: int
    Counts: int
    PSignalAndScNormRatiosAverage: Decimal

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
