from typing import Dict, List, IO


## Считывает файл в словарь, где ключом является заголовок столбца
# (получается из первой строки файла), а значением — список значений в столбце
# @param tableFilename Название файла, из которого происходит чтение
# @param unsafeFlag Разрешить считывание строк, содержащих больше столбцов,
#                   чем есть заголовков
# @param sep Разделитель столбцов
# @returns Словарь, в котором ключами выступают заголовки столбцов, а
#          значениями — список значений в столбце
def ReadTable(tableFilename: str,
              unsafeFlag: bool = False,
              sep: str = '\t') -> Dict[str, List[str]]:
    with open(tableFilename) as inFile:
        return ReadTableFromFileObj(inFile, unsafeFlag, sep)


## Считывает файловый поток в словарь, где ключом является заголовок столбца
# (получается из первой строки файла), а значением — список значений в столбце
# @param inFile Поток для чтения
# @param unsafeFlag Разрешить считывание строк, содержащих больше столбцов,
#                   чем есть заголовков
# @param sep Разделитель столбцов
# @returns Словарь, в котором ключами выступают заголовки столбцов, а
#          значениями — список значений в столбце
def ReadTableFromFileObj(inFile: IO[str],
                         unsafeFlag: bool,
                         sep: str) -> Dict[str, List[str]]:
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
