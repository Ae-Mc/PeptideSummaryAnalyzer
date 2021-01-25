from typing import Dict, List
from .Accession import Accession
from .Sequence import Sequence
from .PeptideTables import PeptideTables
from .PeptideTable import PeptideTable
from decimal import Decimal


class AccessionTables(dict):
    """Хранит список Accession для каждой таблицы и предоставляет функции
    работы с ним

    Словарь, хранящий Accession'ы, распределённые по таблицам.
    Имеет вид: {
        "Номер таблицы": {
            "Имя Accession": Accession
        }
    }

    Args:
        sortedTableNums: отсортированный список номеров таблиц
    """
    sortedTableNums: List[str]

    def __init__(self,
                 seqDB: Dict[str, Sequence],
                 peptideTables: PeptideTables) -> None:
        """См. GetAccessionsPerTable
        """
        self.GetAccessionsPerTable(seqDB, peptideTables)

    def GetAccessionsPerTable(self,
                              seqDB: Dict[str, Sequence],
                              peptideTables: PeptideTables) -> None:
        """Конвертирует PeptideTables в AccessionTables и подсчитывает
        нормализованные значения для каждого Accession

        Получает unused, суммы значений Sc, Precursor Signal и сумму длинн
        последовательностей для каждого Accession, а также нормализованные
        значения Precursor Signal и Sc для каждого файла.
        Этот метод не заполняет поля ScNormToFileNormRatio,
        PSignalNormToFileNormRatio и PSignalAndScNormRatiosAverage!!!

        Args:
            seqDB: словарь с последовательностями, считанными из БД, вида: {
                    "Имя Accession": Sequence
                }
            peptideTables: класс PeptideTables, который будет конвертирован в
                AccessionTables
        """

        self.clear()
        for tableNum, peptideTable in peptideTables.items():
            self[tableNum] = self._GetAccessionsFromTable(peptideTable)
            self._CalculateNormParamsForAccessions(self[tableNum], seqDB)

    def _GetAccessionsFromTable(
            self,
            peptideTable: PeptideTable) -> Dict[str, Accession]:
        """Конвертирует PeptideTables в AccessionTables

        Получает unused, суммы значений Sc, Precursor Signal и сумму длинн
        последовательностей и подсчитывает количество строк с Accession для
        каждого Accession

        Args:
            seqDB: словарь с последовательностями, считанными из БД, вида: {
                    "Имя Accession": Sequence
                }
            peptideTables: класс PeptideTables, который будет конвертирован в
                AccessionTables

        Returns:
            Словарь с Accession'ами вида: {
                "Имя Accession": Accession
            }
        """
        accessions: Dict[str, Accession] = {}
        i = 0
        while i < len(peptideTable):
            curAccession = peptideTable[i].name
            if curAccession not in accessions:
                accessions[curAccession] = Accession(name=curAccession)
            accessions[curAccession].Counts += 1
            accessions[curAccession].Unused = Decimal(peptideTable[i].unused)
            accessions[curAccession].ScSumm += Decimal(peptideTable[i].sc)
            accessions[curAccession].PSignalSumm += (
                Decimal(peptideTable[i].precursorSignal) if
                peptideTable[i].precursorSignal != '' else Decimal(0))
            accessions[curAccession].SeqlenSumm += len(
                peptideTable[i].sequence)
            i += 1
        return accessions

    def _CalculateNormParamsForAccessions(
            self,
            accessions: Dict[str, Accession],
            seqDB: Dict[str, Sequence]) -> None:
        """Подсчитывает нормализованные параметры Sc и Psignal для всех
        Accession'ов из словаря accessions

        Этот метод не заполняет поля ScNormToFileNormRatio,
        PSignalNormToFileNormRatio и PSignalAndScNormRatiosAverage!!!

        Args:
            accessions: словарь вида: {
                    "Имя Accession": Accession
                }
            seqDB: словарь с последовательностями, считанными из БД, вида: {
                    "Имя Accession": Sequence
                }
        """
        for accession, curAccession in accessions.items():
            curAccession.ScNorm = curAccession.ScSumm / seqDB[accession].len
            curAccession.PSignalNorm = (
                curAccession.PSignalSumm / seqDB[accession].len)

    def GenerateAccessionsBunchOverAllTables(
            self) -> Dict[str, Dict[str, Accession]]:
        """Получает словарь с Accession, каждый из которых разбит по таблицам

        Если конкретного Accession нет в конкретной таблице, то эта таблица
        отсутствует в словаре таблиц этого Accession

        Returns:
            Словарь вида: {
                "Имя Accession": {
                    "Номер таблицы": Accession
                }
            }
        """
        accessions: Dict[str, Dict[str, Accession]] = {}
        for tableName, table in self.items():
            for accessionName, accession in table.items():
                if accessionName not in accessions:
                    accessions[accessionName] = {}
                accessions[accessionName][tableName] = accession
        return accessions

    def RemoveAccessionFromAllTables(self, accession: str) -> None:
        """Удаляет Accession из всех таблиц

        Args:
            accession: имя Accession
        """
        for table in self.values():
            table.pop(accession, None)
