process MVP_ANNO {
    container "https://depot.galaxyproject.org/singularity/python:3.9"
    maxRetries 2
    maxForks 3

    input:
    path(invcf)
    path(mvp_tab)

    output:
    path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}" + ".mvp.vcf"
    """
#!/usr/bin/env python

import csv
import re

infile = "${mvp_tab}"
invcf = "${invcf}"
outvcf = '${outvcf}'

mvpdict = dict()

with open(infile, mode = 'rt', newline = '') as incon:
    inreader = csv.reader(incon, delimiter = '\\t')
    for row in inreader:
        if row[0].startswith('#'):
            continue
        chrom = row[0]
        pos = int(row[1])
        ref = row[2]
        alt = row[3]
        aaref = row[4]
        aaalt = row[5]
        gene = row[6]
        txpt = row[7].replace(';', '/')
        score = row[8]
        if chrom not in mvpdict.keys():
            mvpdict[chrom] = dict()
        if pos not in mvpdict[chrom].keys():
            mvpdict[chrom][pos] = dict()
        if alt not in mvpdict[chrom][pos].keys():
            mvpdict[chrom][pos][alt] = list()
        mvpdict[chrom][pos][alt].append((ref, aaref, aaalt, gene, txpt, score))

revcomp = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

try:
    incon = open(invcf, 'rt')
    outcon = open(outvcf, 'wt')
    for line in incon:
        # File header
        if line.startswith('##'):
            outcon.write(line)
            continue
        if line.startswith('#CHROM'):
            outcon.write('##INFO=<ID=MVP,Number=.,Type=String,Description="Gene, transcript, and pathogenicity score from MVP">\\n')
            outcon.write(line)
            continue
        # Data lines
        row = line.split('\\t', 8)
        chrom = row[0].replace('chr', '')
        pos = int(row[1])
        ref = row[3]
        alt = row[4]
        info = row[7]
        # Find the variant
        if chrom not in mvpdict.keys() or pos not in mvpdict[chrom].keys() or len(ref) != 1 or len(alt) != 1 or alt == '*':
            mvpstring = '.'
        else:
            # sanity check to make sure the reference matches, try reverse complement if not
            temp = mvpdict[chrom][pos][list(mvpdict[chrom][pos].keys())[0]]
            assert len(temp) > 0
            ref1 = temp[0][0]
            if ref1 != ref:
                ref = revcomp[ref]
                alt = revcomp[alt]
                if ref1 != ref:
                    print("Mismatched reference allele at chromosome " + chrom + " position " + str(pos))
                    mvpstring = '.'
            elif alt not in mvpdict[chrom][pos].keys():
                mvpstring = '.'
            else:
                scores = mvpdict[chrom][pos][alt]
                assert len(scores) > 0
                # sanity check to make sure the amino acid changes match
                bcsq = re.sub(";.*", "", re.sub(".*;BCSQ=", "", info))
                aaref = [scr[1] for scr in scores]
                aaalt = [scr[2] for scr in scores]
                aa_patterns = ["|\\d+" + aaref[i] + ">\\d+" + aaalt[i] + "|" for i in range(len(scores))]
                scores1 = list()
                for i in range(len(scores)):
                    if re.search(aa_patterns[i], bcsq) != None:
                        scores1.append(scores[i])
                    else:
                        print("Mismatched amino acid change at chromosome " + chrom + " position " + str(pos) + " gene " + scores[i][3])
                if len(scores1) == 0:
                    mvpstring = '.'
                else:
                    mvpstring = ','.join(['|'.join(scr[3:]) for scr in scores])
        # write VCF line with info
        newinfo = info + ';MVP=' + mvpstring
        outcon.write('\\t'.join(row[:7]) + '\\t' + newinfo + '\\t' + row[8])
finally:
    incon.close()
    outcon.close()
    """
}