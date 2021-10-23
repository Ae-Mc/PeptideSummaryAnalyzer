"""См. абстрактный класс ColumnNames."""

from dataclasses import dataclass, fields
from typing import Sequence


@dataclass
class ColumnNames:
    """Хранит имена столбцов и их позиции"""

    def test_column_names(self, column_names: Sequence[str]):
        """Проверяет column_names на соответствие именам столбцов и их
        позициям, прописанным в этом классе

        Args:
            column_names (Sequence[str]): список заголовков столбцов,
                подлежащих проверке

        Raises:
            IndexError: выбрасывается при условии, что количество столбцов
                меньше чем позиция одного из столбцов, прописанных в этом
                классе
            ValueError: выбрасывается если один из столбцов не совпадает
        """
        for field in fields(self):
            position = self.__dict__[field.name][0]
            if position > len(column_names):
                raise IndexError(
                    "Number of columns is smaller than position"
                    f" of column {self.__dict__[field.name][1]}"
                )
            current_column = column_names[position]
            if current_column != self.__dict__[field.name][1]:
                raise ValueError(
                    f"Unknown column on position {position}!"
                    f" Expected {self.__dict__[field.name][1]}"
                    f", got {current_column}"
                )
