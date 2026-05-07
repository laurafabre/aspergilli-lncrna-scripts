import sys
import re

with open(sys.argv[1], "r") as gtf, open(sys.argv[2], "r") as non_coding_list, open(sys.argv[3], "w") as processed_gtf:

    # Read the list of accepted transcript IDs (e.g., XLOC_...)
    non_coding_ids = {}
    for line in non_coding_list:
        line = line.strip()
        xloc_match = re.search(r"(XLOC_[0-9]{6})", line)
        if xloc_match:
            xloc = xloc_match.group(1)
            non_coding_ids[xloc] = line

    # Store the transcript IDs that we want to retain
    kept_transcript_ids = set()

    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue

        feature_type = fields[2]
        attributes = fields[8]

        # Extract transcript_id and class_code
        tid_match = re.search(r'transcript_id "([^"]+)"', attributes)
        class_code_match = re.search(r'class_code "([=ux])"', attributes)

        if feature_type == "transcript":
            if class_code_match:
                class_code = class_code_match.group(1)
                if tid_match:
                    transcript_id = tid_match.group(1)
                    kept_transcript_ids.add(transcript_id)
                    processed_gtf.write(line)
        elif feature_type == "exon":
            if tid_match:
                transcript_id = tid_match.group(1)
                if transcript_id in kept_transcript_ids:
                    processed_gtf.write(line)
