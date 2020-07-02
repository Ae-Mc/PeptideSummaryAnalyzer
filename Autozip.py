<<<<<<< HEAD
#!python
=======
#!/bin/env python3
>>>>>>> InDevelopment
import zipfile
from os import listdir
from os.path import join

with zipfile.ZipFile("PeptideSummaryAnalyzer.zip",
                     'w',
                     zipfile.ZIP_STORED) as archive:
    archive.write("PeptideSummaryAnalyzer.py")
    for filename in listdir("Classes"):
        if not filename.startswith("__") and not filename.startswith("."):
            archive.write(join("Classes", filename))
