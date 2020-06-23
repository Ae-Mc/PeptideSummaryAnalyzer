from typing import Dict, List, Union, Optional
from Classes.Comparable import Comparable
from Classes.Sequence import Sequence


class Input:
    inputPath: str
    seqDB: Dict[str, Sequence]
    unused: Comparable
    contrib: Comparable
    __confPeptide: Comparable
    __confID: Comparable
    __isProteinGroupFilter: bool
    isConfID: bool
    whiteList: Optional[List[str]]
    blackList: Optional[List[str]]
    minGroupsWithAccession: int
    maxGroupAbsence: int
    proteinPilotVersion: str

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
