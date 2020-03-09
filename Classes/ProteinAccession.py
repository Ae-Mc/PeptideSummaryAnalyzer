class ProteinAccession:

    """ Класс для хранения данных об Accession из файлов ProteinSummary """

    name: str
    occurences: int
    unused: float

    def __init__(self, name: str, unused: float, occurences: int = 0):
        self.name = name
        self.unused = unused
        self.occurences = occurences

    def __eq__(self, other):
        if not isinstance(other, ProteinAccession):
            return NotImplemented
        return self.__dict__ == other.__dict__
