"""См. класс ProteinColumns."""

from dataclasses import dataclass
from typing import Tuple
from classes.base_classes import ColumnNames


@dataclass
class ProteinColumns(ColumnNames):

    """Хранит номера столбцов для Protein таблиц
    Attributes:
        n: номер столбца N и его имя
        accession: номер столбца accession и его имя
        unused: номер столбца unused и его имя
    """

    n: Tuple[int, str] = (0, "N")  # pylint: disable=invalid-name
    accession: Tuple[int, str] = (6, "Accession")
    unused: Tuple[int, str] = (1, "Unused")
