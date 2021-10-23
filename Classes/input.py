from enum import Enum
from typing import List, Optional, Tuple
from Classes.comparable import Comparable
from Classes.SequenceDatabase import SequenceDatabase


class ProteinConfidenceType(Enum):
    NONE = 0
    DEFAULT = 1
    VALUE = 2


class Input:
    """Содержит все нужные параметры для запуска обработки

    Attributes:
        root_path: путь для поиска базы данных fasta
        input_path: путь для таблиц
        output_path: путь для выходных файлов
        fdr: параметры FDR фильтра
        seq_db: база данных с длинами последовательностей для Accession
        protein_confidence: параметр фильтра Protein confidence
        protein_grouping_confidence: параметр фильтра Protein grouping (conf)
        conf_peptide: параметр фильтра Peptide confidence
        exclusion_list: чёрный список Accession
        min_groups_with_accession: минимум групп с Accession
        max_group_absence: максимальное количество таблиц в группе без
            Accession
        TODO: rewrite attributes docs
    """

    root_path: str
    inputPath: str
    output_path: str = "Output"
    __fdr: str
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
        return self.__fdr

    def set_fdr(self, value: str):
        value = value.strip()
        if any([value == "default", value == ""]):
            self.__fdr = "default"
        else:
            raise NotImplementedError(
                f"FDR param {value} support not implemented yet"
            )

    @property
    def protein_confidence(self) -> Comparable:
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
        if val.strip() == "default":
            self.protein_confidence_type = ProteinConfidenceType.DEFAULT
        elif len(val.strip()):
            self.protein_confidence_type = ProteinConfidenceType.VALUE
            self.__protein_confidence = Comparable(val)
        else:
            self.protein_confidence_type = ProteinConfidenceType.NONE

    @property
    def protein_grouping_confidence(self) -> Comparable:
        return self.__protein_grouping_confidence

    def set_protein_grouping_confidence(self, val: str):
        self.__protein_grouping_confidence = Comparable(val)

    @property
    def conf_peptide(self) -> Comparable:
        return self.__conf_peptide

    def set_conf_peptide(self, val: str):
        self.__conf_peptide = Comparable(val)
