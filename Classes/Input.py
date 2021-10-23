from Classes.comparable import Comparable
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
    __fdr: str
    seqDB: SequenceDatabase
    __proteinConfidence: Comparable
    __proteinGroupingConfidence: Comparable
    __confPeptide: Comparable
    isProteinConfidence: bool
    blackList: Optional[Tuple[str, List[str]]]
    minGroupsWithAccession: Optional[int]
    maxGroupLack: Optional[int]

    def __init__(self) -> None:
        self.minGroupsWithAccession = None
        self.maxGroupLack = None

    def getFDRStr(self) -> str:
        return self.__fdr

    def setFDR(self, value: str):
        value = value.strip()
        if any([value == "default", value == ""]):
            self.__fdr = "default"
        else:
            raise NotImplementedError(f"FDR param {value} support not implemented yet")

    @property
    def proteinConfidence(self) -> Optional[Comparable]:
        if self.isProteinConfidence:
            return self.__proteinConfidence
        else:
            return None

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
