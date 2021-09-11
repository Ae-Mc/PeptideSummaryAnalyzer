from sqlite3.dbapi2 import Cursor
from typing import List, Tuple


class ParamFilters:
    """Отвечает за фильтрацию по параметрам contribution и confidence

    Attributes:
        cursor: cursor, через который создаются таблицы"""

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    # TODO
    def confidenceDefault(
        self,
        originalTableName: str = "peptide_row",
        searchTableName: str = """peptide_row JOIN peptide_accession
        AS {pseudonym}
        ON peptide_row.id = row_id""",
    ) -> List[Tuple]:
        return self.cursor.execute(
            f"""SELECT DISTINCT table_number, accession FROM (
            SELECT table_number, accession, confidence
            FROM {searchTableName.format(pseudonym="t1")}
            WHERE confidence < 99
            GROUP BY table_number, accession
            HAVING (SELECT COUNT(*)
                FROM {searchTableName.format(pseudonym="t2")}
                WHERE confidence >= 95
                    AND table_number = table_number
                    AND accession = accession
                ) < 2
            );"""
        ).fetchall()
