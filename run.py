from sys import argv
import os

# Arguments:
# N, time, dt, name(optional, determines if filewrite) (optional plotting args)

# Solve WITHOUT writing data:
if len(argv) <= 4:
    # Run simulation:
    os.system("./main.x %s %s %s 0" % (argv[1], argv[2], argv[3]))

# Solve and write data:
else:
    # Check if path present, else make path:
    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists("data/%s" % argv[4]):
        os.mkdir("data/%s" % argv[4])

    # Run simulation:
    os.system("./main.x %s %s %s %s" % (argv[1], argv[2], argv[3], argv[4]))

    # Optional plotting:
    if len(argv) >= 6:
        s = "python plot.py %s" % argv[4]
        
        for i in range(5, len(argv)): # argv[5] and up is plotting commands
            s += " %s" % argv[i]

    os.system(s) # run plot.py

