import re
from sys import argv
fh = open(argv[1])
this_peak = []
print("Gene\tCategory\tCleavage Position\tP-Value\tFragment Abundance\tWeighted Fragment Abundance\tNormalised Weighted Fragment Abundance\tAlignment Score\tShort Read ID\tShort Read Abundance\tNormalised Short Read Abundance\tDuplex")
for line in fh:
    if line != "-----\n":
        line = line.strip("\n")
        this_peak.append(line)
    else:
        this_peak = "\n".join(this_peak)
        csv_peak = this_peak.replace('\t\n', '\t"')
        csv_peak += "\""
        print csv_peak
        this_peak = []

