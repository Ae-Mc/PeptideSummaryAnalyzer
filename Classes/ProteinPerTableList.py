from typing import List

from Classes.ProteinGroup import ProteinGroup
from Classes.ProteinGroupsDB import ProteinGroupsDB


class ProteinPerTableList(dict):
    def __init__(self, db: ProteinGroupsDB) -> None:
        self.ReadFromProteinGroupsDB(db)

    def ReadFromProteinGroupsDB(self, db: ProteinGroupsDB) -> None:
        table: List[ProteinGroup]
        for tableNum, table in db.items():
            self[tableNum] = []
            for group in table:
                if group.representativeAccession not in self[tableNum]:
                    self[tableNum].append(group.representativeAccession)
                for accession in group.accessions:
                    if accession not in self[tableNum]:
                        self[tableNum].append(accession)
