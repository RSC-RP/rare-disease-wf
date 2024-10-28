// Sort VCF (or other file) names based on chromosome names extracted from a FAIDX.
process FILEORDER {
    container "https://depot.galaxyproject.org/singularity/python:3.9"
    label 'process_single'

    input: 
    path(files)      // List of files
    val(chromstring) // String of chromosome names output by CHROMNAMES
    
    output: stdout   // Single string. Space-delimited list of VCF file names.

    script:
    """
#!/usr/bin/env python

import re

chromosomes = \"\"\"${chromstring}\"\"\".strip().split(' ')
chromosomes.append('OTHER')

vcflist = "$files".split(' ')

# Find which file names contain which chromosome names
matches1 = [set([i for i in range(len(chromosomes)) if chromosomes[i] in f]) for f in vcflist]
# Error if any VCFs don't match chromosomes
if any([len(m) == 0 for m in matches1]):
    print([vcflist[i] for i in range(len(vcflist)) if len(matches1[i]) == 0])
    print(chromosomes)
    raise Exception("File doesn't match any chromosome names.")

# Find unique matches
matches2 = [None for f in vcflist]
while any([len(m) == 1 for m in matches1]):
    for i in range(len(chromosomes)):
        fileind = [j for j in range(len(vcflist)) if matches1[j] == {i}]
        if len(fileind) == 0:
            continue
        if len(fileind) > 1:
            print([vcflist[j] for j in fileind])
            print(chromosomes[i])
            raise Exception("Multiple files match single chromosome")
        matches2[fileind[0]] = i
        matches1 = [m - {i} for m in matches1]

# Try to match in the other direction
matches3 = [set([j for j in range(len(vcflist)) if i in matches1[j]]) for i in range(len(chromosomes))]
while any([len(m) == 1 for m in matches3]):
    for j in range(len(vcflist)):
        chromind = [i for i in range(len(chromosomes)) if matches3[i] == {j}]
        if len(chromind) == 0:
            continue
        # If multiple chromosomes match a single file, try to get the one with the longest name.
        if len(chromind) > 1:
            maxlen = max([len(chromosomes[i]) for i in chromind])
            chromind = [i for i in chromind if len(chromosomes[i]) == maxlen]
        if len(chromind) > 1:
            print(vcflist[j])
            print([chromosomes[i] for i in chromind])
            raise Exception("Multiple chromosomes match single file")
        matches2[j] = chromind[0]
        matches3 = [m - {j} for m in matches3]

if all([m == None for m in matches2]):
    print(vcflist)
    print(chromosomes)
    raise Exception("Unable to match chromosomes. May need to improve algorithm.")

# If any are still ambiguous, use file prefix and suffix to sort it out
if any([m == None for m in matches2]):
    nchar_prefix = set()
    ncharp1_suffix = set()
    for j in [j for j in range(len(vcflist)) if matches2[j] != None]:
        chrom = chromosomes[matches2[j]]
        thisfile = vcflist[j]
        matchpos = [(m.start(), m.end()) for m in re.finditer(chrom, thisfile)]
        if len(matchpos) > 1:
            continue
        nchar_prefix.add(matchpos[0][0])
        ncharp1_suffix.add(len(thisfile) - matchpos[0][1])
    if len(nchar_prefix) != 1 or len(ncharp1_suffix) != 1:
        print(vcflist)
        print(chromosomes)
        raise Exception("Unable to match chromosomes. May need to improve algorithm.")
    nchar_prefix = nchar_prefix.pop()
    ncharp1_suffix = ncharp1_suffix.pop()
    for j in [j for j in range(len(vcflist)) if matches2[j] == None]:
        chrom = vcflist[j][nchar_prefix:-ncharp1_suffix]
        try:
            matches2[j] = chromosomes.index(chrom)
        except ValueError:
            print(chromosomes)
            print(vcflist[j])
            raise Exception("Unable to match file to chromosomes.")

# Sort vcf files by index list
vcflist = [f for _, f in sorted(zip(matches2, vcflist))]

print(' '.join(vcflist))
    """
}
