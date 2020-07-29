from typing import Union


def IsFloat(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


# Класс для хранения параметров фильтра
class Comparable:
    op: str
    __val: str
    __isFileToCompare: bool

    def __init__(self, op: str = "", val: str = None) -> None:
        self.__isFileToCompare: bool = True
        if val is not None:
            self.op: str = op
            self.val: str = val
        elif op is not None:
            self.GetComparable(op)
        else:
            self.op = None
            self.val = None

    @property
    def val(self):
        return self.__val

    @val.setter
    def val(self, value):
        self.__val = value
        try:
            float(value)
            self.__isFileToCompare = False
        except (TypeError, ValueError):
            self.__isFileToCompare = True

    def GetComparable(self, paramString: str):
        """ Получение param

        Получение param, общего для всех фалов, либо для каждого своего, из
        строки.
        Формат: [операция][[имя файла] или [число]]
        Примеры:
            >=99
            < paramList.txt"""

        paramString = paramString.strip()
        self.op = ''.join([ch for ch in paramString if ch in "!=<>"])

        paramString = paramString[len(self.op):].strip()
        if IsFloat(paramString):
            self.val = float(paramString)
        elif len(self.op) and len(paramString):
            with open(paramString) as paramStringFile:
                strings = paramStringFile.read().replace(' ', '\t').split('\n')
                paramStrings = {}
                for string in strings:
                    string = string.strip()
                    if len(string):
                        paramStrings[string.split('\t')[0]] = (
                            string.split('\t')[1])
                self.val = paramStrings
        else:
            self.val = None

    def compare(self,
                value: Union[str, float, int],
                filename: str = None) -> bool:
        if self.val is None:
            return True
        elif not self.__isFileToCompare:
            return eval(f"{value}{self.op}{self.val}")
        elif self.val and filename in self.val:
            return eval(f"{value}{self.op}{self.val[filename]}")
        else:
            return True

    def __str__(self):
        return f"{self.op} {self.val}"

    def __repr__(self):
        return f"{self.op} {self.val}"
