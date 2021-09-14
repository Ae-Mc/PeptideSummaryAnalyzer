from pprint import pprint
from sqlite3 import Connection, Cursor
from typing import Any, Iterable, List, Tuple

from Classes.DB.db import DB
from Classes.Functions import FindFastaFile, GetFileLines
from Classes.Input import Input
from Classes.PeptideColumns import PeptideColumns
from Classes.RawPeptideTables import RawPeptideTables
from Classes.SequenceDatabase import SequenceDatabase


class tempDB:

    connection: Connection
    cursor: Cursor

    def getDefaultConfidenceFiltredAccessions(self) -> List[Tuple[str, str]]:
        return self.execute(
            """SELECT DISTINCT table_number, accession FROM (
            SELECT table_number, accession
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = row_id
            WHERE confidence >= 95 AND confidence < 99
            GROUP BY table_number, accession
            HAVING COUNT(*) >= 2
            UNION SELECT table_number, accession
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = row_id
            WHERE confidence >= 99
            GROUP BY table_number, accession);"""
        ).fetchall()

    # TODO
    def addPrecursorAndScoreAccessionSums(self) -> None:
        pass

    # TODO
    def getTableSums(self) -> None:
        self.execute(
            """SELECT
                SUM(precursor_signal) AS precursor_signal_sum,
                SUM(score) AS score_sum"""
        )

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


if __name__ == "__main__":
    inp = Input()
    inp.rootPath = "./Presets/preset_IDexcl"
    inp.inputPath = inp.rootPath + "/Input"
    inp.seqDB = SequenceDatabase.fromFile(FindFastaFile(inp.rootPath))
    peptideTables = RawPeptideTables(PeptideColumns(), inp.inputPath)
    db = DB()

    # Заполняем таблицы
    db.fillers.fillSequence(inp.seqDB)
    db.fillers.fillPeptide(peptideTables)

    db.functions.testFastaDatabase()
    try:
        db.fillers.fillExclusion(GetFileLines(inp.rootPath + "/IDexcl.txt") or [])
    except FileNotFoundError:
        print("ID exclusion file not found")

    db.fdr.default()
    db.functions.applyExclusion()

    print("Количество появлений каждого accession по всем таблицам")
    pprint(
        db.execute(
            """CREATE TABLE overall_accession_count AS
            SELECT accession, COUNT(*) AS count
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = peptide_accession.row_id
            GROUP BY accession
            ORDER BY accession"""
        ).fetchall()
    )
    print("Количество появлений каждого accession в каждой таблице")
    pprint(
        db.execute(
            """CREATE TABLE accession_count_per_table AS
            SELECT table_number, accession, COUNT(*) AS count
            FROM peptide_row JOIN peptide_accession
                ON peptide_row.id = peptide_accession.row_id
            GROUP BY table_number, accession
            ORDER BY table_number, accession"""
        ).fetchall()
    )

    print(
        "Выбор репрезентативных accession для строк и заполнение таблицы"
        " accession_table"
    )
    # TODO
    # pprint(
    #     db.execute(
    #         """SELECT row_id, accession FROM peptide_accession
    #         WHERE """))

    print("Получение accession, не соответствующих confidence default фильтру")
    pprint(db.paramFilters.confidenceDefault())

    pprint(db.getDB())
    db.close()
