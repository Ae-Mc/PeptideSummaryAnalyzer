from typing import Dict, List, IO


def ReadTable(tableFilename: str,
              unsafeFlag: bool = False,
              sep: str = '\t') -> Dict[str, List[str]]:
    """Считывает файл в словарь, где ключом является заголовок столбца
    (получается из первой строки файла), а значением — список значений
    в столбце

    Args:
        tableFilename Название файла, из которого происходит чтение
        unsafeFlag Разрешить считывание строк, содержащих больше столбцов,
                      чем есть заголовков
        sep Разделитель столбцов

    Returns:
        Словарь, в котором ключами выступают заголовки столбцов, а
        значениями — список значений в столбце"""
    with open(tableFilename) as inFile:
        return ReadTableFromFileObj(inFile, unsafeFlag, sep)


def ReadTableFromFileObj(inFile: IO[str],
                         unsafeFlag: bool,
                         sep: str) -> Dict[str, List[str]]:
    """Считывает файловый поток в словарь, где ключом является заголовок
    столбца (получается из первой строки файла), а значением — список значений
    в столбце

    Args:
        inFile Поток для чтения
        unsafeFlag Разрешить считывание строк, содержащих больше столбцов,
                      чем есть заголовков
        sep Разделитель столбцов

    Returns:
        Словарь, в котором ключами выступают заголовки столбцов, а
        значениями — список значений в столбце"""
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
