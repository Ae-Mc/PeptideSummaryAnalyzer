from typing import List, Union, Optional, Tuple
from Classes.Comparable import Comparable
from Classes.SequenceDatabase import SequenceDatabase


class Input:
    """Содержит все нужные параметры для запуска обработки

    Attributes:
        rootPath: путь для поиска базы данных fasta
        inputPath: путь для таблиц
        outputPath: путь для выходных файлов
        fdr: параметры FDR фильтра
        seqDB: база данных с длинами последовательностей для Accession
        proteinConfidence: параметр фильтра Protein confidence
        proteinGroupingConfidence: параметр фильтра Protein grouping (conf)
        confPeptide: параметр фильтра Peptide confidence
        blackList: чёрный список Accession
        minGroupsWithAccession: минимум групп с Accession
        maxGroupAbsence: максимальное количество таблиц в группе без Accession
    """

    rootPath: str
    inputPath: str
    outputPath: str = "Output"
    fdr: str
    seqDB: SequenceDatabase
    __proteinConfidence: Comparable
    __proteinGroupingConfidence: Comparable
    __confPeptide: Comparable
    isProteinConfidence: bool
    blackList: Optional[Tuple[str, List[str]]]
    minGroupsWithAccession: int
    maxGroupAbsence: int

    @property
    def proteinConfidence(self):
        return self.__proteinConfidence

    @proteinConfidence.setter
    def proteinConfidence(self, val):
        self.isProteinConfidence = True
        if val.strip() == "default":
            self.__proteinConfidence = Comparable()
        elif len(val.strip()):
            self.__proteinConfidence = Comparable(val)
        else:
            self.isProteinConfidence = False

    @property
    def proteinGroupingConfidence(self):
        return self.__proteinGroupingConfidence

    @proteinGroupingConfidence.setter
    def proteinGroupingConfidence(self, val):
        self.__proteinGroupingConfidence = Comparable(val)

    @property
    def confPeptide(self):
        return self.__confPeptide

    @confPeptide.setter
    def confPeptide(self, val):
        self.__confPeptide = Comparable(val)
