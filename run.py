from sys import argv
from os import system

# Usage:
# time, dt, name (determine if filewrite), "no" (optional, turns off plotting)

if len(argv) >= 4: # solve and write data
    system("mkdir data/%s" % argv[3])
    system("./main.x %s %s %s" % (argv[1], argv[2], argv[3]))

    flag = True
    if len(argv) >= 5:
        if argv[4] == "no":
            flag = False

    if flag == True: # also automatically plot results
        system("python plot.py %s &" % argv[3])

else: # solve without writing data
    system("./main.x %s %s 0" % (argv[1], argv[2]))
