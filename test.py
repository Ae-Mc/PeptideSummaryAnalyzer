#!python3.7
import os
os.system('python3.7 PeptideSummaryAnalyzer_10.py "" "IDexcl.txt" ' +
          '"y" "EFRA_cont.fasta" ">=0.3" "" "default" 1 3')
ERRORFLAG = False

for filename in os.listdir("./Output"):
    with open(f"OutputOriginal/{filename}") as originalFile,\
         open(f"Output/{filename}") as newFile:
        originalFileContent = originalFile.read()
        newFileContent = newFile.read()
        for i in range(0, len(originalFileContent)):
            if(i == len(newFileContent) or
               originalFileContent[i] != newFileContent[i]):
                print(f"File named {filename} not equal!")
                ERRORFLAG = True
                break

print("Test completed {}!".format("unsuccessful" if ERRORFLAG
                                  else "successful"))
