from sqlite3 import Cursor


class OutputGrouping:
    """Отвечает за работу #Output filter -> groups.

    Attributes:
        cursor: экземляр класса Cursor для связи с базой данных
    """

    cursor: Cursor

    def __init__(self, cursor: Cursor) -> None:
        self.cursor = cursor

    def applyGroupFilter(
        self,
        maxGroupLack: int,
        minGroupsWithAccession: int,
    ) -> None:
        """Применение фильтра по группам к peptide_with_sum

        Таблицы, начинающиеся с одинакового числа входят в одну группу, например:
        1.1, 1.2, 1.3 входят в одну группу, 2.1, 2.2 входят в другую и т. д.
        Из всех таблиц удаляются все Accession, присутствующие в меньше, чем в
        minGroupsWithAccession группах. Accession считается отсутствуюющим в
        группе, если он отсутствует в больше, чем maxGroupAbsence таблицах внутри
        группы.

        Args:
            maxGroupLack: максимально возможное количество таблиц без accession
                в группе
            minGroupsWithAccession: минимальное количество групп с accession
        """

        self.cursor.execute(
            f"""--sql
            WITH accession_bunch AS (
                SELECT DISTINCT *
                FROM (
                    SELECT DISTINCT GET_GROUP_NUMBER(table_number) AS group_number
                    FROM peptide_with_sum
                ) JOIN (SELECT DISTINCT accession FROM peptide_with_sum)
            ),
            table_count_per_group AS (
                SELECT GET_GROUP_NUMBER(table_number) AS group_number, COUNT(*) AS count
                FROM (SELECT DISTINCT table_number FROM peptide_with_sum)
                GROUP BY group_number
            ),
            accession_lack AS (
                SELECT
                    group_number,
                    accession_bunch.accession,
                    (SELECT count FROM table_count_per_group counts
                    WHERE counts.group_number = accession_bunch.group_number)
                    - COUNT(CASE WHEN table_number IS NOT NULL THEN 1 END)
                        AS lack
                FROM accession_bunch LEFT JOIN peptide_with_sum
                    ON group_number = GET_GROUP_NUMBER(table_number)
                        AND accession_bunch.accession = peptide_with_sum.accession
                GROUP BY group_number, accession_bunch.accession
            )
            DELETE FROM peptide_with_sum WHERE accession IN (
                SELECT accession
                FROM accession_lack
                GROUP BY accession
                HAVING COUNT(CASE WHEN lack <= (?) THEN 1 END) < (?)
            );""",
            [maxGroupLack, minGroupsWithAccession],
        )
