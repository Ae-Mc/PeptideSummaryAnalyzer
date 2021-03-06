#!/usr/bin/env python3
import os

if "Output" in os.listdir():
    for filename in os.listdir("./Output"):
        os.remove(f"./Output/{filename}")

os.system('/usr/bin/env python3 PeptideSummaryAnalyzer.py "5" "" '
          '"IDexcl.txt" "y" "EFRA_cont.fasta" ">=0.3" "" "default" 1 3')
ERRORFLAG = False

for filename in os.listdir("./OutputOriginal"):
    with open(f"OutputOriginal/{filename}") as originalFile,\
         open(f"Output/{filename}") as newFile:
        originalFileContent = originalFile.read()
        newFileContent = newFile.read()
        line = 1
        character = 0
        for i in range(0, len(originalFileContent)):
            if newFileContent[i] == '\n':
                line += 1
                character = 0
            character += 1
            if(i == len(newFileContent) or
               originalFileContent[i] != newFileContent[i]):
                print(f"File named {filename} not equal on line {line}" +
                      f" on character {character}!")
                ERRORFLAG = True
                break

print("Test completed {}successful!".format("un" if ERRORFLAG else ""))
