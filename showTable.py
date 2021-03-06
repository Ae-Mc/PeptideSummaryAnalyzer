#!/usr/bin/env python3
from typing import List, Tuple
from sys import argv
from Classes.ProteinTable import ProteinTable
from Classes.BaseClasses.Table import Table as BaseTable


class Table(BaseTable):
    def Load(self, a):
        super().Load(a)


class Column:
    width: int
    height: int
    data: List[str]
    header: str

    def __init__(self,
                 columnName: str = None,
                 columnData: List[str] = None) -> None:
        if columnName is None or columnData is None:
            self.width = 0
            self.height = 0
        else:
            self.Load(columnName, columnData)

    def Load(self,
             columnName: str,
             columnData: List[str]) -> None:
        self.header = columnName
        self.data = columnData
        self.height = len(columnData)
        self.width = len(columnName)
        for value in columnData:
            if len(value) > self.width:
                self.width = len(value)


class TablePrinter:
    columns: List[Column]
    columnPadding: int
    topLeftCorner: str
    topRightCorner: str
    downLeftCorner: str
    downRightCorner: str

    topT: str
    downT: str
    leftT: str
    rightT: str

    verticalSymbol: str
    horizontalSymbol: str

    crossSymbol: str

    def __init__(self,
                 table: BaseTable = None,
                 filename: str = None,
                 topLeftCorner: str = '┌',
                 topRightCorner: str = '┐',
                 downLeftCorner: str = '└',
                 downRightCorner: str = '┘',
                 topT: str = '┬',
                 downT: str = '┴',
                 leftT: str = '├',
                 rightT: str = '┤',
                 verticalSymbol: str = '│',
                 horizontalSymbol: str = '—',
                 crossSymbol: str = '┼'):
        self.columns: List[Column] = []
        self.columnPadding = 1
        self.topLeftCorner: str = topLeftCorner
        self.topRightCorner: str = topRightCorner
        self.downLeftCorner: str = downLeftCorner
        self.downRightCorner: str = downRightCorner

        self.topT: str = topT
        self.downT: str = downT
        self.leftT: str = leftT
        self.rightT: str = rightT

        self.verticalSymbol: str = verticalSymbol
        self.horizontalSymbol: str = horizontalSymbol

        self.crossSymbol: str = crossSymbol

        if table is not None:
            self.LoadTable(table)
        elif filename is not None:
            self.LoadFile(filename)

    def LoadFile(self, filename: str):
        table = Table(filename, unsafeFlag=True)
        self.LoadTable(table)

    def LoadTable(self, table: BaseTable):
        for i, column in enumerate(table[0]):
            self.columns.append(Column(
                column, [table[val][i] for val in range(1, len(table))]))

    def Print(self, columnPadding: int = None):
        if columnPadding is not None:
            self.SetColumnPadding(columnPadding)

        # Print headers string
        self.PrintDelimiterLine(columnDelimiterSymbol=self.topT,
                                leftSymbol=self.topLeftCorner,
                                rightSymbol=self.topRightCorner)
        print(end=self.verticalSymbol)
        for column in self.columns:
            print(f"{' '*self.columnPadding}" +
                  f"{column.header.ljust(column.width + self.columnPadding)}",
                  end=self.verticalSymbol)
        print()
        self.PrintDelimiterLine()
        for i in range(0, self.columns[0].height):
            print(end=self.verticalSymbol)
            for column in self.columns:
                print(f"{' '*self.columnPadding}" +
                      "{}".format(
                          column.data[i].ljust(column.width +
                                               self.columnPadding)),
                      end=self.verticalSymbol)
            print()
        self.PrintDelimiterLine(columnDelimiterSymbol=self.downT,
                                leftSymbol=self.downLeftCorner,
                                rightSymbol=self.downRightCorner)

    def PrintDelimiterLine(self,
                           columnDelimiterSymbol: str = None,
                           leftSymbol: str = None,
                           rightSymbol: str = None):
        if columnDelimiterSymbol is None:
            columnDelimiterSymbol = self.crossSymbol
        if leftSymbol is None:
            leftSymbol = self.leftT
        if rightSymbol is None:
            rightSymbol = self.rightT
        print(end=leftSymbol)
        for column in self.columns[:-1]:
            print('─' * (column.width + self.columnPadding * 2),
                  end=columnDelimiterSymbol)
        print('─' * (self.columns[-1].width + self.columnPadding * 2),
              end=rightSymbol + '\n')

    def SetColumnPadding(self, columnPadding: int):
        self.columnPadding = columnPadding


if len(argv) > 1:
    try:
        columnPadding = int(argv[1])
        filenames = argv[2:]
    except ValueError:
        columnPadding = 1
        filenames = argv[1:]
    for filename in filenames:
        columnNums: Tuple[int, ...] = tuple()
        table = Table(filename, unsafeFlag=True)
        if filename.endswith("ProteinSummary.txt"):
            columnNums = (0, 1, 6, 8, 9,)
        elif filename.endswith("PeptideSummary.txt"):
            columnNums = (0, 1, 6, 10, 11, 12, 22)
        if len(columnNums):
            for i, line in enumerate(table):
                table[i] = [line[j] for j in columnNums]
        tablePrinter = TablePrinter(table=table)
        tablePrinter.Print(columnPadding=columnPadding)
