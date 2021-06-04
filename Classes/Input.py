from typing import List, Union, Optional, Tuple
from Classes.Comparable import Comparable
from Classes.SequenceDatabase import SequenceDatabase


class Input:
    """Содержит все нужные параметры для запуска обработки

    Attributes:
        rootPath: путь для поиска базы данных fasta
        inputPath: путь для таблиц
        outputPath: путь для выходных файлов
        seqDB: база данных с длинами последовательностей для Accession
        confID: параметр фильтра confID
        confPeptide: параметр фильтра confPeptide
        blackList: чёрный список Accession
        minGroupsWithAccession: минимум групп с Accession
        maxGroupAbsence: максимальное количество таблиц в группе без Accession
    """

    rootPath: str
    inputPath: str
    outputPath: str = "Output"
    seqDB: SequenceDatabase
    __confPeptide: Comparable
    __proteinConfidence: Comparable
    __isProteinGroupFilter: bool
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
    def confPeptide(self):
        return self.__confPeptide

    @confPeptide.setter
    def confPeptide(self, val):
        self.__confPeptide = Comparable(val)

    @property
    def isProteinGroupFilter(self):
        return self.__isProteinGroupFilter

    @isProteinGroupFilter.setter
    def isProteinGroupFilter(self, value: Union[str, bool]):
        if isinstance(value, str):
            value = value.lower()
            self.__isProteinGroupFilter = True if value == "y" else False
        else:
            self.__isProteinGroupFilter = value
