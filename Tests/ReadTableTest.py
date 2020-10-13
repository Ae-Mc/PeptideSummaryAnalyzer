from Classes.ReadTable import ReadTableFromFileObj
from io import StringIO
import unittest


class ReadTableTest(unittest.TestCase):

    def testSafeReadTableException(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\n" +
            "5\t6\t7\t8\t9\n" +
            "a\tb\tc\td\n")
        with self.assertRaises(IndexError):
            ReadTableFromFileObj(inputFile, False, '\t')

    def testUnsafeReadTable(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\t\n" +
            "5\t6\t7\t8\t\n" +
            "a\tb\tc\td\n")
        self.assertDictEqual(
            ReadTableFromFileObj(inputFile, True, '\t'),
            {
                'A': ["1", "5", "a"],
                'B': ["2", "6", "b"],
                'C': ["3", "7", "c"],
                'D': ["4", "8", "d"]
            })

    def testSafeReadTable(self):
        inputFile = StringIO(
            "A\tB\tC\tD\n" +
            "1\t2\t3\t4\n" +
            "5\t6\t7\t8\n" +
            "a\tb\tc\td\n")
        self.assertDictEqual(
            ReadTableFromFileObj(inputFile, False, '\t'),
            {
                'A': ["1", "5", "a"],
                'B': ["2", "6", "b"],
                'C': ["3", "7", "c"],
                'D': ["4", "8", "d"]
            })
