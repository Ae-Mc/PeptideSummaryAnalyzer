from sqlite3 import Cursor


class Summaries:
    """Отвечает за заполнение таблицы peptide_with_sum

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def fillPeptideWithSum(self) -> list:
        return self.cursor.execute(
            """--sql
            WITH norm_sums_per_table AS (
                SELECT
                    table_number,
                    SUM(sc_norm) AS sc_norm_sum,
                    SUM(peptide_intensity_norm) AS peptide_intensity_norm_sum
                FROM (
                    SELECT
                        table_number,
                        SUM(score) / LENGTH(sequence.sequence) AS sc_norm,
                        SUM(peptide_intensity) / LENGTH(sequence.sequence)
                            AS peptide_intensity_norm
                    FROM peptide JOIN sequence
                        ON sequence.accession = peptide.accession
                    GROUP BY table_number, peptide.accession
                )
                GROUP BY table_number
            )
            INSERT INTO peptide_with_sum (
                table_number,
                accession,
                count,
                sc_sum,
                peptide_intensity_sum,
                sc_norm,
                peptide_intensity_norm,
                sc_norm_to_file_norm_ratio,
                peptide_intensity_norm_to_file_norm_ratio
            )
            SELECT
                table_number,
                peptide.accession AS accession,
                COUNT(*) AS count,
                SUM(score) AS sc_sum,
                SUM(peptide_intensity) AS peptide_intensity_sum,
                SUM(score) / LENGTH(sequence.sequence) AS sc_norm,
                SUM(peptide_intensity) / LENGTH(sequence.sequence)
                    AS peptide_intensity_norm,
                SUM(score)
                    / LENGTH(sequence.sequence)
                    / (SELECT sc_norm_sum FROM norm_sums_per_table
                        WHERE norm_sums_per_table.table_number = peptide.table_number)
                    AS sc_norm_to_file_norm_ratio,
                SUM(peptide_intensity)
                    / LENGTH(sequence.sequence)
                    / (SELECT peptide_intensity_norm_sum FROM norm_sums_per_table
                        WHERE norm_sums_per_table.table_number = peptide.table_number)
                    AS peptide_intensity_norm_to_file_norm_ratio
            FROM peptide JOIN sequence ON sequence.accession = peptide.accession
            GROUP BY table_number, peptide.accession;"""
        ).fetchall()
