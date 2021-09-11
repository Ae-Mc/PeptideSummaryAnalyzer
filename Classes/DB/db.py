from Classes.DB.Fillers import Fillers
from Classes.DB.FDR import FDR
from Classes.DB.ParamFilters import ParamFilters
from typing import Any, Iterable, List, Tuple
from Classes.DB.Creators import Creators
from Classes.Functions import IsReversed
from sqlite3.dbapi2 import Connection, Cursor, connect


class DB:
    """Отвечает за создание базы данных и функций.

    Содержит другие классы дл взаимодействия с БД.

    Attributes:
        connection: соединение с БД
        cursor: connection.cuesor()
        initializers: класс отвечающий за создание таблиц
        paramFilters: класс, отвечающий за фильтрацию по параметрам confidence
            и contribution
        fdr: класс с FDR фильтрами"""

    connection: Connection
    cursor: Cursor
    initializers: Creators
    paramFilters: ParamFilters
    fdr: FDR

    def __init__(self) -> None:
        self.connection = connect(":memory:")
        self.cursor = self.connection.cursor()
        self.connection.create_function("IS_REVERSED", 1, IsReversed)
        self.initializers = Creators(self.cursor)
        self.initializers.createAllTables()
        self.fillers = Fillers(self)
        self.paramFilters = ParamFilters(self.cursor)
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

    def close(self) -> None:
        self.connection.close()
