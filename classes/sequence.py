"""См. класс Sequence"""

from dataclasses import dataclass


@dataclass
class Sequence:
    """Хранит последовательность

    Attributes:
        accession: имя accession, к которому относится данная
            последовательность
        len: длина последовательности
        desc: описание последовательности
        seq: обработанная последовательность, которая включает только буквы
        raw_seq: сама последовательность
    """

    accession: str = ""
    desc: str = ""
    __raw_seq: str = ""
    __seq: str = ""

    def __init__(self, accession: str = "", desc: str = "", rawSeq: str = ""):
        self.accession = accession
        self.desc = desc
        self.raw_seq = rawSeq

    @property
    def seq(self) -> str:
        """Обработанная последовательность"""
        return self.__seq

    @property
    def len(self) -> int:
        """Длина обработанной последовательности"""
        return len(self.seq)

    @property
    def raw_seq(self) -> str:
        """Сама последовательность"""
        return self.__raw_seq

    @raw_seq.setter
    def raw_seq(self, value: str):
        self.__seq = "".join(
            [ch for ch in value if ch in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
        )
        self.__raw_seq = value

    def __repr__(self):
        return f"{self.seq} ({self.len})"
