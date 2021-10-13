from dataclasses import dataclass


@dataclass
class Sequence:
    """Хранит последовательность

    Attributes:
        accession: имя accession, к которому относится данная последовательность
        len: длина последовательности
        desc: описание последовательности
        seq: обработанная последовательность, которая включает только буквы
        rawSeq: сама последовательность
    """

    accession: str = ""
    desc: str = ""
    __rawSeq: str = ""
    __seq: str = ""

    def __init__(self, accession: str = "", desc: str = "", rawSeq: str = ""):
        self.accession = accession
        self.desc = desc
        self.rawSeq = rawSeq

    @property
    def seq(self) -> str:
        return self.__seq

    @property
    def len(self) -> int:
        return len(self.seq)

    @property
    def rawSeq(self) -> str:
        return self.__rawSeq

    @rawSeq.setter
    def rawSeq(self, value: str):
        self.__seq = "".join([ch for ch in value if ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"])
        self.__rawSeq = value

    @property
    def __repr__(self):
        return f"{self.seq} ({self.len})"
