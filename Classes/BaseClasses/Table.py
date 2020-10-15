from abc import ABC, abstractmethod
from typing import List, Any


class Table(ABC, list):
    """Читает таблицу и возвращает её как список строк, где каждая строка —
    кортеж значений строки

    Attributes:
        unsafeFlag: разрешает строкам таблицы быть разной длины
    """
    unsafeFlag: bool

    def __init__(self,
                 tableFilename: str = None,
                 unsafeFlag: bool = False):
        self.unsafeFlag = unsafeFlag
        if tableFilename:
            self.Load(tableFilename)

    @abstractmethod
    def Load(self,
             tableFilename: str) -> List[Any]:
        self.clear()
        with open(tableFilename) as inFile:
            lines = inFile.read().split('\n')
        if len(lines):
            baseLen = len(lines[0].split("\t"))
            for line in filter(lambda x: len(x), lines):
                values = tuple(line.split("\t"))
                if not self.unsafeFlag and len(values) != baseLen:
                    raise ValueError("Table with variable len lines detected!"
                                     f" File: {tableFilename}")
                self.append(values)
        return self
