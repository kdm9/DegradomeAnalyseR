# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the
# GPLv3 License
from sys import argv

fh = open(argv[1])
this_peak = []
columns = [
        "Gene",
        "Category",
        "Cleavage Position",
        "P-Value",
        "Fragment Abundance",
        "Weighted Fragment Abundance",
        "Normalised Weighted Fragment Abundance",
        "Alignment Score",
        "Short Read ID",
        "Short Read Abundance",
        "Normalised Short Read Abundance",
        "Duplex"
        ]

print "\t".join(columns)

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
