process FILTER_PATHOGENIC {
    container "https://depot.galaxyproject.org/singularity/python:3.9"

    // Default filtering thresholds
    //params.mvp_thresh = 0.9
    //params.cadd_thresh = 7
    //params.dann_thresh = 0.999
    //params.gerp_thresh = 5.5

    input:
    path(invcf)

    output:
    path(outvcf)

    script:
    outvcf = "${invcf.simpleName}" + ".pathogenic.vcf"

"""
#!/usr/bin/env python

try:
    incon = open("$invcf", mode = 'rt')
    outcon = open("$outvcf", mode = 'wt')

    for line in incon:
        if line.startswith('#'):
            outcon.write(line)
            continue
        fields = line.split('\\t', 8)
        info = fields[7].split(';')
        if 'MetaSVM_pred=D' in info or 'FATHMM_pred=D' in info or \
        'Polyphen2_HDIV_pred=D' in info or 'Polyphen2_HVAR_pred=D' in info or \
        'fathmm-MKL_coding_pred=D' in info or \
        'ExonicFunc.refGene=frameshift_insertion' in info or \
        'ExonicFunc.refGene=frameshift_deletion' in info or \
        'Func.refGene=splicing' in info:
            outcon.write(line)
            continue

        mvp_fields = [x for x in info if x.startswith('MVP=')]
        if len(mvp_fields) == 1 and mvp_fields[0] != 'MVP=.':
            mvp_items = mvp_fields[0].split(',')
            scores = [float(x.split('|')[-1]) for x in mvp_items]
            if max(scores) > $params.mvp_thresh:
                outcon.write(line)
                continue

        cadd_fields = [x for x in info if x.startswith('CADD_raw=')]
        if len(cadd_fields) == 1 and cadd_fields[0] != 'CADD_raw=.' and \
        float(cadd_fields[0][9:]) > $params.cadd_thresh:
            outcon.write(line)
            continue

        dann_fields = [x for x in info if x.startswith('DANN_score=')]
        if len(dann_fields) == 1 and dann_fields[0] != 'DANN_score=.' and \
        float(dann_fields[0][11:]) >= $params.dann_thresh:
            outcon.write(line)
            continue

        gerp_fields = [x for x in info if x.startswith('GERP++_RS=')]
        if len(gerp_fields) == 1 and gerp_fields[0] != 'GERP++_RS=.' and \
        float(gerp_fields[0][10:]) > $params.gerp_thresh:
            outcon.write(line)
            continue
        
        cln_fields = [x for x in info if x.startswith('CLNSIG=')]
        if len(cln_fields) == 1 and cln_fields[0] != 'CLNSIG=.' and \
        'athogenic' in cln_fields[0]:
            outcon.write(line)
            continue
        
finally:
    incon.close()
    outcon.close()
"""
}
