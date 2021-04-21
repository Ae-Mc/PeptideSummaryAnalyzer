from dataclasses import dataclass, fields
from typing import Sequence


@dataclass
class ColumnNames:

    def TestColumnNames(self, columnNames: Sequence[str]):
        for field in fields(self):
            position = self.__dict__[field.name][0]
            if position > len(columnNames):
                raise IndexError("Number of columns is smaller than position"
                                 f" of column {self.__dict__[field.name][1]}")
            curColumn = columnNames[position]
            if curColumn != self.__dict__[field.name][1]:
                raise ValueError(f"Unknown column on position {position}!"
                                 f' Expected "{self.__dict__[field.name][1]}"'
                                 f', got "{curColumn}"')
