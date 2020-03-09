from Classes.ProteinAccession import ProteinAccession


class ProteinAccessionWithMaxUnused(ProteinAccession):
    # Если максимальный unused у Accession сразу в нескольких группах
    ambiguity: bool

    def __init__(
            self,
            name: str,
            unused: float,
            occurences: int = 0,
            ambiguity: bool = False) -> None:
        self.name = name
        self.unused = unused
        self.occurences = occurences
        self.ambiguity = ambiguity

    def __eq__(self, other):
        if not isinstance(other, ProteinAccessionWithMaxUnused):
            return NotImplemented
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return (f"{self.name} (U:{self.unused}, Occurs:{self.occurences}, "
                f"Ambiguity:{self.ambiguity})")
