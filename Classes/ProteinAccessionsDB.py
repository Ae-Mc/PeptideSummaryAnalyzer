from typing import Dict, List
from decimal import Decimal
from .Sequence import Sequence
from .ProteinTable import ProteinTable
from .ProteinAccession import ProteinAccession
from .BaseClasses.ProteinDB import ProteinDB


class ProteinAccessionsDB(ProteinDB):
    """Хранит базу данных с Protein Accession в формате: {
        "номер таблицы": {
            "имя Accession": ProteinAccession
        }
    }

    Вычисляет максимальные параметры Unused для каждого Accession.
    Служит для поиска репрезентативного Accession среди списка Accession.
    """

    def LoadFromTables(self, proteinTables: Dict[str, ProteinTable]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц

        Args:
            proteinTables: словарь вида: {
                    "номер таблицы": ProteinTable
                }
        """
        for tableNum, table in proteinTables.items():
            curUnused: Decimal = Decimal(0)
            i = 0
            while i < len(table):
                if table[i].unused != Decimal(0):
                    curUnused = table[i].unused
                accession = ProteinAccession(table[i].name, curUnused)
                if accession.name in self:
                    self[accession.name].unused = max(
                        self[accession.name].unused, accession.unused)
                    self[accession.name].occurences += 1
                else:
                    self[accession.name] = accession
                i += 1

    def GetRepresentative(self,
                          accessions: List[str],
                          seqDB: Dict[str, Sequence]) -> str:
        """Получает репрезентативный Accession для списка Accession

        Args:
            accessions: список имён Accession
            seqDB: база данных последовательностей Accession

        Returns:
            Имя репрезентативного Accession
        """
        representativeAccession = ProteinAccession("", Decimal(0), 0)
        for accession in sorted(accessions):
            if(
                self[accession].unused > representativeAccession.unused or
                (self[accession].unused == representativeAccession.unused and
                 self[accession].occurences >
                 representativeAccession.occurences) or
                (self[accession].unused == representativeAccession.unused and
                 self[accession].occurences ==
                 representativeAccession.occurences and
                 seqDB[accession].len >
                 seqDB[representativeAccession.name].len)):
                representativeAccession = self[accession]
        return representativeAccession.name
