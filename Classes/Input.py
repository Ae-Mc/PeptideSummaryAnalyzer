from typing import List, Union, Optional
from Classes.Comparable import Comparable
from Classes.SequenceDatabase import SequenceDatabase


class Input:
    """Содержит все нужные параметры для запуска обработки

    Attributes:
        inputPath: путь для таблиц
        outputPath: путь для выходных файлов
        seqDB: база данных с длинами последовательностей для Accession
        unused: параметр фильтра unused
        contrib: параметр фильтра contrib
        confID: параметр фильтра confID
        confPeptide: параметр фильтра confPeptide
        whiteList: белый список Accession
        blackList: чёрный список Accession
        minGroupsWithAccession: минимум групп с Accession
        maxGroupAbsence: максимальное количество таблиц в группе без Accession
        proteinPilotVersion: версия ProteinPilot (от этого зависит формат
            таблиц)
    """
    inputPath: str
    outputPath: str = "Output"
    seqDB: SequenceDatabase
    unused: Comparable
    contrib: Comparable
    __confPeptide: Comparable
    __confID: Comparable
    __isProteinGroupFilter: bool
    isConfID: bool
    blackList: Optional[List[str]]
    minGroupsWithAccession: int
    maxGroupAbsence: int

    @property
    def confID(self):
        return self.__confID

    @confID.setter
    def confID(self, val):
        self.isConfID = True
        if val.strip() == "default":
            self.__confID = Comparable(None, None)
        elif len(val.strip()):
            self.__confID = Comparable(val)
        else:
            self.isConfID = False

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
            self.__isProteinGroupFilter = True if value == 'y' else False
        else:
            self.__isProteinGroupFilter = value
