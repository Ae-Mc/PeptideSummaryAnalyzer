from typing import List

from Classes.ProteinGroup import ProteinGroup
from Classes.ProteinGroupsDB import ProteinGroupsDB


class ProteinPerTableList(dict):
    """Словарь для хранения Protein таблиц в виде: {
        "номер таблицы": ["Accession1", "Accession2", ..., "AccessionN"]
    }
    """
    def __init__(self, db: ProteinGroupsDB) -> None:
        """См. ReadFromProteinGroupsDB"""
        self.ReadFromProteinGroupsDB(db)

    def ReadFromProteinGroupsDB(self, db: ProteinGroupsDB) -> None:
        """Считывает ProteinGroupsDB в словарь

        Словарь вида: {
            "номер таблицы": ["Accession1", "Accession2", ..., "AccessionN"]
        }

        Args:
            db: ProteinGroupsDB из которой происходит считывание
        """
        table: List[ProteinGroup]
        for tableNum, table in db.items():
            self[tableNum] = []
            for group in table:
                if group.representativeAccession not in self[tableNum]:
                    self[tableNum].append(group.representativeAccession)
                for accession in group.accessions:
                    if accession not in self[tableNum]:
                        self[tableNum].append(accession)
