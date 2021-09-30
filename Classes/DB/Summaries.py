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
            -- INSERT INTO peptide_with_sum (
            --     table_number,
            --     accession,
            --     count,
            --     sc_sum,
            --     peptide_intensity_sum,
            --     sc_norm,
            --     peptide_intensity_norm,
            --     sc_norm_to_file_norm_ratio,
            --     peptide_intensity_norm_to_file_norm_ratio
            -- )
            SELECT
                table_number,
                peptide.accession,
                COUNT(*),
                SUM(score) AS sc_sum,
                SUM(peptide_intensity) AS peptide_intensity_sum,
                SUM(score) / LENGTH(sequence.sequence) AS sc_norm,
                SUM(peptide_intensity) / LENGTH(sequence.sequence)
                    AS peptide_intensity_norm,
                SUM(score)
                    / LENGTH(sequence.sequence)
                    / (SELECT SUM(sc_norm) FROM (
                            SELECT
                                table_number,
                                SUM(score) / LENGTH(sequence.sequence) AS sc_norm
                            FROM peptide JOIN sequence
                                ON sequence.accession = peptide.accession
                            GROUP BY table_number, peptide.accession
                        ) t1
                        WHERE t1.table_number = peptide.table_number
                        GROUP BY table_number
                    )
                    AS sc_norm_to_file_norm_ratio,
                SUM(peptide_intensity)
                    / LENGTH(sequence.sequence)
                    / (SELECT SUM(peptide_intensity_norm) FROM (
                            SELECT
                                table_number,
                                SUM(peptide_intensity) / LENGTH(sequence.sequence)
                                    AS peptide_intensity_norm
                            FROM peptide JOIN sequence
                                ON sequence.accession = peptide.accession
                            GROUP BY table_number, peptide.accession
                        ) t1
                        WHERE t1.table_number = peptide.table_number
                        GROUP BY table_number
                    )
                    AS peptide_intensity_norm_to_file_norm_ratio
            FROM peptide JOIN sequence ON sequence.accession = peptide.accession
            GROUP BY table_number, peptide.accession;"""
        ).fetchall()

    def createFileSumTable(self) -> None:
        pass
