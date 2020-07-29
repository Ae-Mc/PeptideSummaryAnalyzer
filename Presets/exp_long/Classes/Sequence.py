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

    def __repr__(self):
        return f"{self.seq} ({self.len})"
