from sys import argv
import os

# Arguments:
# time, dt, name(optional, determines if filewrite) (optional plotting args)

# Solve WITHOUT writing data:
if len(argv) <= 3:
    os.system("./main.x %s %s 0" % (argv[1], argv[2])) # run simulation

# Solve and write data:
else:
    # Check if path present, else make path:
    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists("data/%s" % argv[3]):
        os.mkdir("data/%s" % argv[3])
        
    os.system("./main.x %s %s %s" % (argv[1], argv[2], argv[3]))# run simulation

    # Optional plotting:
    if len(argv) >= 5:
        s = "python plot.py %s" % argv[3]
        
        for i in range(4, len(argv)): # argv[4] and up is plotting commands
            s += " %s" % argv[i]

    os.system(s) # run plot.py

