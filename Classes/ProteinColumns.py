from typing import Tuple
from .BaseClasses.ColumnNames import ColumnNames, dataclass


@dataclass
class ProteinColumns(ColumnNames):

    """Хранит номера столбцов для Protein таблиц

    Attributes:
        accession: номер столбца accession и его имя
        unused: номер столбца unused и его имя
    """

    accession: Tuple[int, str] = (6, "Accession")
    unused: Tuple[int, str] = (1, "Unused")
