from typing import Union


class Comparable:
    """Класс для хранения параметров фильтра

    Позволяет хранить параметры фильтра и сравнивать их с любыми значениями

    Attributes:
        op: операция сравнения
        val: значение, с которым будет производится сравнение, либо словарь
            вида: {
                "имя файла": значение
            }.
            Также может быть равно None, в этом случае все операции сравнения
            вернут истину
    """
    op: str
    __val: str
    __isFileToCompare: bool

    def __init__(self, op: str = "", val: str = None) -> None:
        self.__isFileToCompare: bool = True
        if val is not None:
            self.op = op
            self.val = val
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
        try:
            float(value)
            self.__val = value
            self.__isFileToCompare = False
        except (TypeError, ValueError):
            if (value is None
                    or (isinstance(value, str) and len(value.strip()) == 0)):
                self.__val = value
            elif isinstance(value, dict):
                self.__val = value
                self.__isFileToCompare = True
            else:
                self._LoadFromFile(value)
                self.__isFileToCompare = True

    def GetComparable(self, paramString: str):
        """Получение параметра

        Получение параметра, общего для всех фалов, либо для каждого своего, из
        строки.
        Формат: [операция][[имя файла] или [число]]

        Args:
            paramString: строка, которая будет ковертирована в параметр

        Примеры:
            >=99
            < paramList.txt"""

        paramString = paramString.strip()
        self.op = ''.join([ch for ch in paramString if ch in "!=<>"])

        paramValue = paramString[len(self.op):].strip()
        self.val = paramValue

    def _LoadFromFile(self, filename: str):
        """Получение значений параметра из файла

        Args:
            filename: имя файла, из которого будут считываться значения
        """
        with open(filename) as paramStringFile:
            strings = paramStringFile.read().replace(' ', '\t').split('\n')
            paramStrings = {}
            for string in strings:
                string = string.strip()
                if len(string):
                    paramStrings[string.split('\t')[0]] = (
                        string.split('\t')[1])
            self.val = paramStrings

    def compare(self,
                value: Union[str, float, int],
                filename: str = None) -> bool:
        """Сравнивает значение с параметром

        Args:
            value: значение, которое будет сравниваться
            filename: имя файла, если параметр был считан из файла
        Returns:
            Результат сравнения значения с параметром
        """
        if self.val is None:
            return True
        elif not self.__isFileToCompare:
            return eval(f"{value}{self.op}{self.val}")
        elif self.val and filename in self.val:
            return eval(f"{value}{self.op}{self.val[filename]}")
        return True

    def __str__(self):
        return f"{self.op} {self.val}"

    def __repr__(self):
        return f"{self.op} {self.val}"
