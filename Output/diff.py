import os
for filename in os.listdir():
    if(not filename.split('.')[0].endswith("Original")
       and filename.endswith(".txt")):
        with open(filename) as newFile:
            with open(filename.split('.')[0] + "Original.txt") as origFile:
                newFileContent = newFile.read()
                origFileContent = origFile.read()
                diff = newFileContent.find('\n') - origFileContent.find('\n')
                for i in range(newFileContent.find('\n'), len(newFileContent)):
                    if i - diff + 1 > len(origFileContent):
                        print("{} left in file {}".format(newFileContent[i:],
                                                          filename))
                        break

                    if newFileContent[i] != origFileContent[i - diff]:
                        print("{}!={} on position {} in file {}".format(
                            newFileContent[i], origFileContent[i - diff],
                            i, filename))
                        break
input("Press Enter to exit")
