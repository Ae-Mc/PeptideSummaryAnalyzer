"""Отвечает за работу Protein Grouping фильтра"""
from sqlite3 import Cursor

from Classes.comparable import Comparable


class ProteinGrouping:
    """Отвечает за работу Protein Grouping

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def create_filtred_accession_table(self, confidence: Comparable) -> None:
        """Создаёт промежуточную таблицу с применённым фильтром по confidence.

        Если фильтр confidence не задан (confidence.op is None), то всё равно
        создаёт таблицу, но без применения фильтра."""

        if confidence.operation is not None:
            self.cursor.execute(
                f"""--sql
                CREATE TABLE filtered_peptide_accession AS SELECT *
                FROM peptide_joint
                WHERE confidence {confidence.operation} {confidence.val};"""
            )
        else:
            self.cursor.execute(
                """--sql
                CREATE TABLE filtered_peptide_accession AS
                    SELECT * FROM peptide_joint;"""
            )

    def fill_count_table(self):
        """Заполняет таблицу с количеством появлений каждого accession в каждой
        таблице."""

        self.cursor.execute(
            """--sql
            INSERT INTO accession_count_per_table (
                table_number, accession, count
            )
            SELECT table_number, accession, COUNT(*) AS count_per_table
            FROM filtered_peptide_accession
            GROUP BY table_number, accession
            ORDER BY table_number, accession;
            """
        )

    def create_general_count_table(self):
        """Создаёт таблицу с количеством файлов с accession для каждого
        accession."""

        self.cursor.execute(
            """--sql
            CREATE TABLE general_accession_count AS
            SELECT accession, COUNT(*) AS general_count
            FROM (
                SELECT DISTINCT table_number, accession
                FROM filtered_peptide_accession
            )
            GROUP BY accession
            ORDER BY accession;"""
        )

    def apply_protein_grouping(self) -> None:
        """Заполняет таблицы representative и accession_group."""

        self._create_groups()
        self.cursor.execute(
            """--sql
            WITH first_criteria_representative AS (
                SELECT DISTINCT representative_id, accession
                FROM (
                    SELECT
                        accession,
                        count,
                        MAX(count)
                            OVER(PARTITION BY table_number, representative_id)
                            AS max_count_in_table_in_group,
                        representative_id
                    FROM accession_group
                        JOIN accession_count_per_table USING(accession)
                )
                WHERE count = max_count_in_table_in_group
            ),
            second_criteria_representative AS (
                SELECT representative_id, accession
                FROM (
                    SELECT
                        representative_id,
                        accession,
                        general_count,
                        MAX(general_count) OVER(PARTITION BY representative_id)
                            AS max_count_in_group
                    FROM first_criteria_representative
                        JOIN general_accession_count USING(accession)
                )
                WHERE general_count = max_count_in_group
            ),
            third_criteria_representative AS (
                SELECT representative_id, accession
                FROM (
                    SELECT
                        representative_id,
                        accession,
                        LENGTH(sequence) AS seqlen,
                        MAX(LENGTH(sequence))
                            OVER(PARTITION BY representative_id)
                            AS max_seqlen_in_group
                    FROM second_criteria_representative
                        JOIN sequence USING(accession)
                )
                WHERE seqlen = max_seqlen_in_group
            ),
            fourth_criteria_representative AS (
                SELECT representative_id, accession
                FROM (
                    SELECT
                        representative_id,
                        accession,
                        MAX(accession) OVER(PARTITION BY representative_id)
                            AS max_accession_in_group
                    FROM second_criteria_representative
                )
                WHERE accession = max_accession_in_group
            )
            UPDATE representative
            SET
                representative = (
                    SELECT accession FROM fourth_criteria_representative
                    WHERE representative.id = representative_id
                );
            """
        )

    def _create_groups(self) -> None:
        """Создаёт группы accession'ов. Заполняет таблицу accession_group и создаёт
        записи в таблице representative, в которых в дальнейшем будет заполнено
        поле representative."""

        accession: str
        for [accession] in self.cursor.execute(
            """--sql
            SELECT DISTINCT accession FROM peptide_accession
            ORDER BY accession;
            """
        ).fetchall():
            if (
                self.cursor.execute(
                    """--sql
                    SELECT COUNT(*) FROM accession_group WHERE accession=(?);
                    """,
                    [accession],
                ).fetchone()[0]
                == 0
            ):
                self.cursor.execute(
                    """INSERT INTO representative DEFAULT VALUES;"""
                )
                representative_id: int = self.cursor.execute(
                    "SELECT last_insert_rowid();"
                ).fetchone()[0]
                self.cursor.execute(
                    """--sql
                    INSERT INTO accession_group (representative_id, accession)
                        SELECT DISTINCT (?), accession FROM peptide_accession
                        WHERE row_id IN (
                            SELECT row_id FROM peptide_accession
                            WHERE accession = (?)
                        )
                        ORDER BY accession
                    ;""",
                    [representative_id, accession],
                )

    def fill_peptide_table(self) -> None:
        """Заполняет таблицу peptide на основе данных из таблиц representative и
        peptide_row."""

        self.cursor.execute(
            """--sql
        INSERT INTO peptide (
            table_number,
            accession,
            confidence,
            score,
            peptide_intensity,
            sequence
        )
        SELECT
            table_number,
            representative,
            confidence,
            score,
            peptide_intensity,
            sequence
        FROM peptide_row JOIN (
            SELECT DISTINCT row_id, representative
            FROM peptide_accession JOIN accession_group USING(accession)
                JOIN representative ON representative.id = representative_id
        ) ON id = row_id;"""
        ).fetchall()
