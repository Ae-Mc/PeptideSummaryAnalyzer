from typing import Dict, List, IO


def ReadTable(tableFilename: str,
              unsafeFlag: bool = False,
              sep: str = '\t') -> Dict[str, List[str]]:
    """ Считываем файл в словарь, где ключом является заголовок, а
    значением — список значений в столбце.

    Если unsafeFlag выставлен в True, то во всех строках после разбиения по
    табуляции элементы выходящие за пределы списка columns (его длина равна
    количеству заголовков) будут удаляться.
    """

    with open(tableFilename) as inFile:
        return ReadTableFromFileObj(inFile, sep, unsafeFlag)


def ReadTableFromFileObj(inFile: IO[str], sep: str, unsafeFlag: bool):
    strings = inFile.read().split('\n')
    inFile.close()
    table: Dict[str, List[str]] = {}
    columns = strings[0].split(sep)

    for column in columns:
        table[column] = []

    for string in strings[1:]:
        if len(string.strip()):
            i = 0
            for value in string.split(sep):
                if unsafeFlag:
                    if i < len(columns):
                        table[columns[i]].append(value)
                    else:
                        break
                else:
                    table[columns[i]].append(value)
                i += 1

    return table
