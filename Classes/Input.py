from typing import Dict, List, Union
from Classes.Comparable import Comparable
from Classes.Sequence import Sequence


class Input:
    inputPath: str
    seqDB: Dict[str, Sequence]
    unused: Comparable
    contrib: Comparable
    __conf: Comparable
    isDefaultConf: bool
    isProteinGroupFilter: bool
    whiteList: Union[List[str], None]
    blackList: Union[List[str], None]
    minGroupsWithAccession: int
    maxGroupAbsence: int
    proteinPilotVersion: str

    @property
    def conf(self):
        return self.__conf

    @conf.setter
    def conf(self, val):
        self.__conf = Comparable(val.strip())
        if val.strip() == "default":
            self.isDefaultConf = True
        else:
            self.isDefaultConf = False
