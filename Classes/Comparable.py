from re import match, sub
from typing import Optional


class Comparable:
    """Класс для хранения параметров фильтра

    Позволяет хранить параметры фильтра и сравнивать их с любыми значениями

    Attributes:
        op: операция сравнения
        val: значение, с которым будет производится сравнение"""

    op: Optional[str]
    val: Optional[str]

    def __init__(self, op: str = None, val: str = None) -> None:
        if not (val is None or op is None):
            self.op = op
            self.val = val
        elif op is not None:
            self.GetComparable(op)
        else:
            self.op = None
            self.val = None

    def GetComparable(self, paramString: str):
        """Получение параметра

        Формат: [операция][число]

        Args:
            paramString: строка, которая будет ковертирована в параметр

        Примеры:
            >=99"""

        x = match(r"^(([<>]=?)|=) *([-\d]+)$", sub("$==", "=", paramString.strip()))
        if x:
            self.op = x.group(1)
            self.val = x.group(3)
        elif paramString.strip() == "":
            self.op = None
            self.val = None
        else:
            raise ValueError("Неверный формат параметра для сравнения")

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if not (self.op is None or self.val is None):
            return f"{self.op} {self.val}"
        else:
            return ""
