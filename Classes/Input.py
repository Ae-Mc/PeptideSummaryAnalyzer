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
        shouldExtractSequences: нужно ли извлекать последовательности
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
    maxGroupLack: int
    shouldExtractSequences: bool

    def __init__(self) -> None:
        self.shouldExtractSequences = True

    @property
    def proteinConfidence(self) -> Comparable:
        return self.__proteinConfidence

    def setProteinConfidence(self, val: str):
        self.isProteinConfidence = True
        if val.strip() == "default":
            self.__proteinConfidence = Comparable()
        elif len(val.strip()):
            self.__proteinConfidence = Comparable(val)
        else:
            self.isProteinConfidence = False

    @property
    def proteinGroupingConfidence(self) -> Comparable:
        return self.__proteinGroupingConfidence

    def setProteinGroupingConfidence(self, val: str):
        self.__proteinGroupingConfidence = Comparable(val)

    @property
    def confPeptide(self) -> Comparable:
        return self.__confPeptide

    def setConfPeptide(self, val: str):
        self.__confPeptide = Comparable(val)
