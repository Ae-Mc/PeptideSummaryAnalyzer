"""Содержит классы, отвечающие за хранение входных параметров скрипта."""

from enum import Enum
from typing import List, Optional, Tuple
from classes.comparable import Comparable
from classes.sequence_database import SequenceDatabase


class FDRtype(Enum):
    """Типы FDR фильтра"""

    NONE = 0
    DEFAULT = 1


class ProteinConfidenceType(Enum):
    """Типы #Protein filter -> Peptide confidence (value or default)."""

    NONE = 0
    DEFAULT = 1
    VALUE = 2


class Input:
    """Содержит все нужные параметры для запуска обработки

    Attributes:
        root_path: путь для поиска базы данных fasta
        input_path: путь для таблиц
        output_path: путь для выходных файлов
        fdr_type: тип fdr фильтра
        seq_db: база данных с длинами последовательностей для Accession
        protein_confidence_type: тип фильтра Protein confidence (см. класс
            ProteinConfidenceType)
        protein_confidence: параметр фильтра Protein confidence
        protein_grouping_confidence: параметр фильтра Protein grouping (conf)
        conf_peptide: параметр фильтра Peptide confidence
        exclusion_list: чёрный список Accession
        min_groups_with_accession: минимум групп с Accession
        max_group_absence: максимальное количество таблиц в группе без
            Accession
    """

    root_path: str
    inputPath: str
    output_path: str = "Output"
    fdr_type: FDRtype
    seq_db: SequenceDatabase
    protein_confidence_type: ProteinConfidenceType
    __protein_confidence: Comparable
    __protein_grouping_confidence: Comparable
    __conf_peptide: Comparable
    exclusion_list: Optional[Tuple[str, List[str]]]
    min_groups_with_accession: Optional[int]
    max_group_lack: Optional[int]

    def __init__(self) -> None:
        self.min_groups_with_accession = None
        self.max_group_lack = None

    def get_fdr_str(self) -> str:
        """Возвращает значение параметра fdr в виде строки."""
        if self.fdr_type == FDRtype.NONE:
            return ""
        if self.fdr_type == FDRtype.DEFAULT:
            return "default"
        raise NotImplementedError("Unknown fdr_type")

    def set_fdr(self, value: str):
        """Устанавливает значение fdr для фильтра в зависимости от входной
        строки.

        Args:
            value (str): входная строка.

        Raises:
            NotImplementedError: вызывается в случае неизвестного формата
                входной строки.
        """
        value = value.strip()
        if value == "":
            self.fdr_type = FDRtype.NONE
        elif value == "default":
            self.fdr_type = FDRtype.DEFAULT
        else:
            raise NotImplementedError(
                f"FDR param {value} support not implemented yet"
            )

    @property
    def protein_confidence(self) -> Comparable:
        """Возвращает значение protein_confidence, если было задано числовое
        значение protein_confidence.


        Raises:
            ValueError: выбрасывает если protein_confidence не было задано
                числовое значение.

        Returns:
            Comparable: protein_confidence
        """
        if self.protein_confidence_type == ProteinConfidenceType.NONE:
            raise ValueError(
                "Can't get protein confidence if protein confidence unset"
            )
        if self.protein_confidence_type == ProteinConfidenceType.DEFAULT:
            raise ValueError(
                "Can't get protein confidence"
                " if protein confidence set to default"
            )
        return self.__protein_confidence

    def set_protein_confidence(self, val: str) -> None:
        """Устанавливает значение protein_confidence_type и protein_confidence
        (опционально) в зависимости от входной строки.

        Args:
            val (str): входная строка.
        """
        if val.strip() == "default":
            self.protein_confidence_type = ProteinConfidenceType.DEFAULT
        elif len(val.strip()):
            self.protein_confidence_type = ProteinConfidenceType.VALUE
            self.__protein_confidence = Comparable(val)
        else:
            self.protein_confidence_type = ProteinConfidenceType.NONE

    @property
    def protein_grouping_confidence(self) -> Comparable:
        """Значение для фильтра Protein Grouping (conf)."""
        return self.__protein_grouping_confidence

    def set_protein_grouping_confidence(self, val: str):
        """Устанавливает значение для фильтра Protein Grouping (conf) в
        зависимости от входной строки.

        Args:
            val (str): входная строка.
        """
        self.__protein_grouping_confidence = Comparable(val)

    @property
    def conf_peptide(self) -> Comparable:
        """Значение для фильтра Peptide confidence (value)."""
        return self.__conf_peptide

    def set_conf_peptide(self, val: str):
        """Устанавливает значение для фильтра Peptide confidence (value) в
        зависимости от входной строки.

        Args:
            val (str): входная строка.
        """
        self.__conf_peptide = Comparable(val)
