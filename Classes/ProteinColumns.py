from dataclasses import dataclass


@dataclass
class ProteinColumns:

    """Хранит номера столбцов для Protein таблиц

    Attributes:
        accession: номер столбца accession
        unused: номер столбца unused
    """

    accession: int = 6
    unused: int = 1
