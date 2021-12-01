"""См. класс Comparable."""

from re import match, sub
from typing import Optional


class Comparable:
    """Класс для хранения параметров фильтра

    Позволяет хранить параметры фильтра и сравнивать их с любыми значениями

    Attributes:
        operation: операция сравнения
        val: значение, с которым будет производится сравнение"""

    operation: Optional[str]
    val: Optional[str]

    def __init__(self, op: str = None, val: str = None) -> None:
        if not (val is None or op is None):
            self.operation = op
            self.val = val
        elif op is not None:
            self._get_comparable(op)
        else:
            self.operation = None
            self.val = None

    def _get_comparable(self, param_string: str):
        """Получение параметра

        Формат: [операция][число]

        Args:
            paramString: строка, которая будет ковертирована в параметр

        Примеры:
            >=99"""

        matched = match(
            r"^(([<>]=?)|=) *(-?\d+(\.\d+)?)$",
            sub("$==", "=", param_string.strip()),
        )
        if matched:
            self.operation = matched.group(1)
            self.val = matched.group(3)
        elif param_string.strip() == "":
            self.operation = None
            self.val = None
        else:
            raise ValueError("Неверный формат параметра для сравнения")

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if not (self.operation is None or self.val is None):
            return f"{self.operation}{self.val}"
        return ""
