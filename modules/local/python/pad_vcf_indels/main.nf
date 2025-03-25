process PAD_VCF_INDELS {
    shell '/bin/bash', '-euo', 'pipefail'
    container 'https://depot.galaxyproject.org/singularity/pyfaidx:0.7.1--pyh5e36f6f_0'

    input:
        path(invcf)
        path(fasta)
        path(fai)
    output: path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}.padded.vcf"
    println "Padding malformed deletions in VCF"
    println "Output file: ${outvcf}"
    """
#!/usr/bin/env python

from pyfaidx import Fasta

vcffilein = "${invcf}"
vcffileout = "${outvcf}"
fasta = "${fasta}"

genome = Fasta(fasta)

try:
    incon = open(vcffilein, 'r')
    outcon = open(vcffileout, 'w')
    for line in incon:
        # comment lines
        if line.startswith('#'):
            outcon.write(line)
            continue
        fields = line.split('\t', 6) # don't bother splitting up genotypes
        alt = fields[4].split(',')
        # variant lines with properly formatted alleles
        if not any([x == '' for x in alt]):
            outcon.write(line)
            continue
        # variant lines with malformed deletions
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        assert genome[chrom][pos - 1].seq.upper() == ref[0].upper(), "VCF and reference don't match"
        # update fields
        pos -= 1
        add = genome[chrom][pos - 1].seq.upper()
        ref = add + ref
        alt = [add + x for x in alt]
        fields[1] = str(pos)
        fields[3] = ref
        fields[4] = ','.join(alt)
        lineout = '\t'.join(fields)
        outcon.write(lineout)
        
finally:
    incon.close()
    outcon.close()
    """
}
