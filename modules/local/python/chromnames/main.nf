// Extract chromosome names from a FAIDX
process CHROMNAMES {
    container "https://depot.galaxyproject.org/singularity/python:3.9"

    input: path(fai)
    
    output: stdout // Single string. Space-delimited list of chromosome names.

    script:
    """
#!/usr/bin/env python

import csv

chroms = []

with open("$fai", mode = 'rt', newline = '') as incon:
    inreader = csv.reader(incon, delimiter = '\\t')
    for row in inreader:
        chroms.append(row[0])

print(' '.join(chroms))
    """
}
