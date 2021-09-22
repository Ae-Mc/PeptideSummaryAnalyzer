from Classes.DB.ProteinGrouping import ProteinGrouping
from Classes.DB.Functions import Functions
from Classes.DB.Fillers import Fillers
from Classes.DB.FDR import FDR
from typing import Any, Iterable, List, Tuple
from Classes.DB.Creators import Creators
from sqlite3.dbapi2 import Connection, Cursor, connect


class DB:
    """Отвечает за создание базы данных и функций.

    Содержит другие классы дл взаимодействия с БД.

    Attributes:
        connection: соединение с БД
        cursor: connection.cuesor()
        initializers: класс отвечающий за создание таблиц
        proteinGrouping: класс, отвечающий за работу protein grouping фильтра
        fdr: класс с FDR фильтрами"""

    connection: Connection
    cursor: Cursor
    initializers: Creators
    functions: Functions
    fillers: Fillers
    paramFilters: ProteinGrouping
    fdr: FDR

    def __init__(self) -> None:
        self.connection = connect(":memory:")
        self.cursor = self.connection.cursor()
        self.initializers = Creators(self.cursor)
        self.initializers.createAllFunctions()
        self.initializers.createAllTables()
        self.functions = Functions(self.cursor)
        self.fillers = Fillers(self.cursor)
        self.proteinGrouping = ProteinGrouping(self.cursor)
        self.fdr = FDR(self.cursor)

    def getDB(self) -> List[Tuple[Any]]:
        return self.execute(
            """SELECT row_id, table_number, accession, confidence
            FROM
                (SELECT id, table_number, confidence
                 FROM peptide_row) AS peptide_row
                JOIN peptide_accession ON row_id = peptide_row.id"""
        ).fetchall()

    def execute(self, sql: str, parameters: Iterable[Any] = []) -> Cursor:
        return self.cursor.execute(sql, parameters)

    def prettyFetchAll(self, sql: str, parameters: Iterable[Any] = []) -> list:
        cursor = self.cursor.execute(sql, parameters)
        result = cursor.fetchall()
        result.insert(0, tuple(map(lambda x: x[0], cursor.description)))
        return result

    def close(self) -> None:
        self.connection.close()
