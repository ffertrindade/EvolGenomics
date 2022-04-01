import os
import sys

list_file = open(sys.argv[1], "r")
tree_file = sys.argv[2]

for i in list_file: 
    i = i.strip()
    print("ete3 evol -t " + tree_file + " --alg " + i + " -i " +  i.split(".")[0] + ".png -o " + i.split(".")[0] + " -v3 --models bsA1 bsA --leaves --resume | tee " +  i.split(".")[0] + ".stdout.txt")