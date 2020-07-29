from decimal import Decimal


class ProteinAccession:

    """ Класс для хранения данных об Accession из файлов ProteinSummary """

    name: str
    occurences: int
    unused: Decimal

    def __init__(self, name: str, unused: Decimal, occurences: int = 1):
        self.name = name
        self.unused = unused
        self.occurences = occurences

    def __eq__(self, other):
        if not isinstance(other, ProteinAccession):
            return NotImplemented
        return self.__dict__ == other.__dict__

    def __str__(self):
        return (f"{self.name} (unused: {self.unused}, "
                f"occurences: {self.occurences})")

    def __repr__(self):
        return (f"{self.name} (unused: {self.unused}, "
                f"occurences: {self.occurences})")
