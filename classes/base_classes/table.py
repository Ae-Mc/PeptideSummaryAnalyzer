"""См. абстрактный класс Table."""

from abc import ABC, abstractmethod
from typing import List, Any


class Table(ABC, list):
    """Читает таблицу и возвращает её как список строк, где каждая строка —
    кортеж значений строки

    Attributes:
        unsafe_flag: разрешает строкам таблицы быть разной длины
    """

    unsafe_flag: bool

    def __init__(self, tableFilename: str = None, unsafe_flag: bool = False):
        self.unsafe_flag = unsafe_flag
        super().__init__()
        if tableFilename is not None:
            self.load(tableFilename)

    @abstractmethod
    def load(self, table_filename: str) -> List[Any]:
        """Метод загрузки таблицы по умолчанию. Предполагается, что он будет
        переопределён дочерними классами.

        Args:
            table_filename (str): имя файла, в котором хранится таблица

        Raises:
            ValueError: выбрасывается в случае, если unsafe_flag == False и
                были обнаружены строки с разным количеством столбцов

        Returns:
            List[Any]: список кортежей из строк таблицы, разделённых на столбцы
        """
        self.clear()
        with open(table_filename, encoding="utf-8") as in_file:
            lines = in_file.read().split("\n")
        if len(lines):
            base_len = len(lines[0].split("\t"))
            for line in filter(len, lines):
                values = tuple(line.split("\t"))
                if not self.unsafe_flag and len(values) != base_len:
                    raise ValueError(
                        "Table with variable len lines detected!"
                        f" File: {table_filename}"
                    )
                self.append(values)
        return self
