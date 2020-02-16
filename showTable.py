from typing import List, Dict
from sys import argv
from Classes.ReadTable import ReadTable


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


class Table:
    columns: List[Column]
    columnPadding: int

    def __init__(self,
                 table: Dict[str, List[str]] = None,
                 filename: str = None):
        self.columns: List[Column] = []
        self.columnPadding = 1
        if table is not None:
            self.LoadTable(table)
        elif filename is not None:
            self.LoadFile(filename)

    def LoadFile(self, filename: str):
        table = ReadTable(filename)
        self.LoadTable(table)

    def LoadTable(self, table: Dict[str, List[str]]):
        for columnHeader, column in table.items():
            self.columns.append(Column(columnHeader, column))

    def Print(self, columnPadding: int = None):
        if columnPadding is not None:
            self.SetColumnPadding(columnPadding)

        # Print headers string
        self.PrintBlankTableLine()
        print(end='|')
        for column in self.columns:
            print(f"{' '*self.columnPadding}" +
                  f"{column.header.ljust(column.width + self.columnPadding)}",
                  end='|')
        print()
        self.PrintBlankTableLine()
        for i in range(0, self.columns[0].height):
            print(end='|')
            for column in self.columns:
                print(f"{' '*self.columnPadding}" +
                      "{}".format(
                          column.data[i].ljust(column.width +
                                               self.columnPadding)),
                      end='|')
            print()
        self.PrintBlankTableLine()

    def SetColumnPadding(self, columnPadding: int):
        self.columnPadding = columnPadding

    def PrintBlankTableLine(self):
        print('+', end='')
        for column in self.columns:
            print('—' * (column.width + self.columnPadding * 2), end='+')
        print()


if len(argv) > 1:
    for filename in argv[1:]:
        if filename.endswith("PeptideSummary.txt"):
            fileTable = ReadTable(filename)
            columnsToDelete = [
                "N", "Total", "%Cov", "%Cov(50)", "%Cov(95)", "Used", "Names",
                "Annotation", "Modifications", "Cleavages", "dMass", "Prec MW",
                "Prec m/z", "Theor MW", "Theor m/z", "Theor z", "Spectrum",
                "Time", "PrecursorElution"]
            for columnName in columnsToDelete:
                del fileTable[columnName]
            for i in range(0, len(fileTable["Accessions"])):
                fileTable["Accessions"][i] = (
                    fileTable["Accessions"][i].split(';')[0])
            table = Table(table=fileTable)
        else:
            table = Table(filename=filename)
        table.Print()
