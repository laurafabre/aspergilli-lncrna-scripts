import sys
import os
import re

with open(sys.argv[1], "r") as gtf, open(sys.argv[2], "r") as non_coding_list, open(sys.argv[3], "w") as processed_gtf:

    id_table = []

    non_coding_ids = {}
    for line in non_coding_list:
        line = line.rstrip()
        xloc = re.findall(r"XLOC_[0-9]{6}", line)
        non_coding_ids[xloc[0]] = line

    for line in gtf:
        if line.startswith("#"):
            continue
        line = line.rstrip().split("\t")
        if line[2] == "transcript":
            col9 = line[8].split(";")
            class_code = None
            for entry in col9:
                key, value = entry.strip().split(" ")
                if key == "class_code":
                    class_code = value.strip('"')
                    break

            if class_code in ("=", "u", "x"):
                processed_gtf.write("%s\n" % ("\t".join(line)))

