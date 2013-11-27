from sys import argv
from os import system,mkdir

# Usage:
# time, dt, name (determine if filewrite), "no" (optional, turns off plotting)

if len(argv) >= 4: # solve and write data
    # Check if path present
    if not os.path.exists("data/%s" % argv[3]):
        # If not, create dir
        mkdir("data/%s" % argv[3])
        
    system("./main.x %s %s %s" % (argv[1], argv[2], argv[3]))

    flag = 1 # default
    if len(argv) >= 5:
        if argv[4] == "no":
            flag = 0
        elif argv[4] == "film" or argv[4] == "movie":
            flag = 2

    if flag == 1: # automatically plot results
        system("python plot.py %s &" % argv[3])
    if flag == 2: # automatically save a movie
        system("python plot.py %s movie &" % argv[3])

else: # solve without writing data
    system("./main.x %s %s 0" % (argv[1], argv[2]))
