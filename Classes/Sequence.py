from dataclasses import dataclass


@dataclass
class Sequence:
    """Хранит последовательность

    Attributes:
        len: длина последовательности
        desc: описание последовательности
        seq: сама последовательность
    """
    accession: str = ""
    len: int = 0
    desc: str = ""
    __seq: str = ""

    def __init__(self, accession: str = "", desc: str = "", seq: str = ""):
        self.accession = accession
        self.desc = desc
        self.seq = seq

    @property
    def seq(self):
        return self.__seq

    @seq.setter
    def seq(self, seq: str):
        self.len = len(seq)
        self.__seq = seq

    def __repr__(self):
        return f"{self.seq} ({self.len})"
