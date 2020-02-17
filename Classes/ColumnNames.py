from typing import Dict


class ColumnNames:

    accession: str = "Accessions"
    precursorSignal: str = "PrecursorSignal"
    sc: str = "Sc"
    unused: str = "Unused"
    sequence: str = "Sequence"
    confidence: str = "Conf"
    contribution: str = "Contrib"

    def __init__(self,
                 columnsDict: Dict[str, str] = None,
                 **kwargs) -> None:
        if columnsDict is not None:
            for columnClassicName, columnName in columnsDict.items():
                self.__dict__[columnClassicName] = columnName

        for columnClassicName, columnName in kwargs.items():
            self.__dict__[columnClassicName] = columnName
