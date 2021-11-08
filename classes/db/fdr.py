"""См. класс FDR."""

from sqlite3.dbapi2 import Cursor


class FDR:
    """Отвечает за работу FDR фильтра.

    Attributes:
        cursor: экземляр класса Cursor для связи с БД.
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def none(self) -> None:
        """Удаляет все перевёрнутые accession."""
        self.cursor.execute(
            """--sql
            DELETE FROM peptide_accession WHERE IS_REVERSED(accession);
            """
        )
        self.cursor.execute(
            """--sql
            DELETE FROM peptide_row WHERE row_id NOT IN (
                SELECT DISTINCT row_id FROM peptide_joint
            );
            """
        )

    def default(self) -> None:
        """Применяет default FDR фильтр."""

        self.cursor.execute(
            """--sql
            DELETE FROM peptide_accession WHERE IS_REVERSED(accession);"""
        )
        for row_id, table_num in self.cursor.execute(
            """--sql
            SELECT peptide_row.row_id, table_number
            FROM peptide_row LEFT JOIN peptide_accession
                ON peptide_row.row_id = peptide_accession.row_id
            WHERE peptide_accession.row_id IS NULL
            ORDER BY peptide_row.row_id;"""
        ).fetchall():
            self.cursor.execute(
                """--sql
                DELETE FROM peptide_row
                WHERE row_id >= (?) AND table_number = (?);
                """,
                [row_id, table_num],
            )

    def fill_main_accession(self) -> None:
        """Заполняет столбец main_accession в таблице protein_group"""

        self.cursor.execute(
            """--sql
            UPDATE protein_group
            SET main_accession = (
                SELECT accession FROM (
                    SELECT
                        group_id,
                        accession,
                        row_number() OVER (
                            PARTITION BY group_id
                            ORDER BY unused DESC, IS_REVERSED(accession) DESC
                        ) AS row_priority
                    FROM protein_row
                ) pr
                WHERE (
                    row_priority = 1
                    AND protein_group.group_id = pr.group_id
                )
            );
            """
        )

    def initialize_fdr_summary_table(self) -> None:
        """Создаёт строки для всех ProteinSummary файлов в таблице fdr_summary.
        В результате все столбцы будут полностью заполнены NULL, кроме столбца
        table_number."""
        self.cursor.execute(
            """
            --sql
            INSERT INTO fdr_summary (table_number)
            SELECT DISTINCT table_number FROM protein_group;
            """
        )

    def fill_target_and_decoy_columns(self) -> None:
        """Заполняет столбцы target_count и decoy_count в fdr_summary."""
        self.cursor.execute(
            """
            --sql
            UPDATE fdr_summary
            SET
                target_count = COALESCE(
                    (
                        SELECT COUNT(
                            CASE IS_REVERSED(accession) WHEN FALSE THEN 1 END
                        )
                        FROM protein_group JOIN protein_row USING(group_id)
                        WHERE (
                            main_accession = accession
                            AND unused BETWEEN 0.05 AND 0.1
                            AND protein_group.table_number
                                = fdr_summary.table_number
                        )
                    ),
                    0
                ),
                decoy_count = COALESCE(
                    (
                        SELECT COUNT(
                            CASE IS_REVERSED(accession) WHEN TRUE THEN 1 END
                        )
                        FROM protein_group JOIN protein_row USING(group_id)
                        WHERE (
                            main_accession = accession
                            AND unused BETWEEN 0.05 AND 0.1
                            AND protein_group.table_number
                                = fdr_summary.table_number
                        )
                    ),
                    0
                );
            """
        )

    def fill_fdr_data(self, k: float = None) -> None:
        """Заполняет таблицу fdr_data"""
        placeholders = (
            ("k", "JOIN fdr_summary USING (table_number)")
            if k is None
            else ("(?)", "")
        )
        self.cursor.execute(
            f"""--sql
            INSERT INTO fdr_data (
                table_number,
                table_n,
                unused,
                accession,
                index_n,
                observed_fdr
            ) SELECT
                table_number,
                table_n,
                unused,
                accession,
                index_n,
                index_n * 100 * {placeholders[0]} / table_n AS observed_fdr
            FROM (
                    SELECT
                        table_number,
                        N AS table_n,
                        unused,
                        accession,
                        row_number() OVER (PARTITION BY table_number)
                            AS index_n
                    FROM protein_group
                        JOIN protein_row USING (group_id)
                    WHERE
                        IS_REVERSED(main_accession) = TRUE
                        AND main_accession = accession
                )
                {placeholders[1]};
            """,
            [] if k is None else [k],
        )

    def fill_a_and_b_columns(self, table_range: int) -> None:
        """Заполняет столбцы a и b в fdr_summary."""

        self.cursor.execute(
            """
            --sql
            WITH temp_params AS (
                SELECT
                    table_number,
                    X,
                    Y,
                    Z,
                    G,
                    (Y * Z - G * X) / (n * Z - X * X) AS T
                FROM (
                    SELECT
                        table_number,
                        SUM(table_n) AS X,
                        SUM(LN(observed_fdr)) AS Y,
                        SUM(table_n * table_n) AS Z,
                        SUM(table_n * LN(observed_fdr)) AS G,
                        COUNT(*) AS n
                    FROM fdr_data
                    WHERE index_n <= (?)
                    GROUP BY table_number
                )
            )
            UPDATE fdr_summary
            SET
                a = (
                    SELECT EXP(T) FROM temp_params
                    WHERE temp_params.table_number = fdr_summary.table_number
                ),
                b = (
                    SELECT (G - X * T) / Z FROM temp_params
                    WHERE temp_params.table_number = fdr_summary.table_number
                );
            """,
            [table_range],
        )

    def fill_additional_statistic_columns(self, table_range: int) -> None:
        """Заполняет столбцы squared_R, MAE и MAPE в таблице fdr_summary."""

        self.cursor.execute(
            """
            --sql
            WITH fdr_params AS (
                SELECT table_number, a * EXP(b * table_n) AS FDRc, observed_fdr
                FROM fdr_data JOIN fdr_summary USING(table_number)
                WHERE index_n <= (?)
            ),
            mO AS (
                SELECT table_number, SUM(observed_fdr) / (?) AS mO
                FROM fdr_data
                WHERE index_n <= (?)
                GROUP BY table_number
            )
            UPDATE fdr_summary
            SET
                squared_R = (
                    SELECT
                        1 - SUM((observed_fdr - FDRc) * (observed_fdr - FDRc))
                        / SUM((observed_fdr - mO) * (observed_fdr - mO))
                    FROM fdr_params JOIN mO USING(table_number)
                    WHERE mO.table_number = fdr_summary.table_number
                ),
                MAE = (
                    SELECT SUM(observed_fdr - FDRc) / (?)
                    FROM fdr_params JOIN mO USING(table_number)
                    WHERE mO.table_number = fdr_summary.table_number
                ),
                MAPE = (
                    SELECT SUM(ABS(observed_fdr - FDRc) / FDRc) / (?)
                    FROM fdr_params JOIN mO USING(table_number)
                    WHERE mO.table_number = fdr_summary.table_number
                )
            ;
            """,
            [table_range, table_range, table_range, table_range, table_range],
        )

    def apply_fdr_k(self, fdr: float) -> None:
        """Применяет FDR k/kest диапазон фильтр."""

        for table_number, row_id in self.cursor.execute(
            """--sql
            SELECT table_number, LN((?) / a) / b FROM fdr_summary;
            """,
            [fdr],
        ).fetchall():
            self.cursor.execute(
                """
                --sql
                DELETE FROM peptide_row
                WHERE N >= (?) AND table_number = (?);
                """,
                [row_id, table_number],
            )
