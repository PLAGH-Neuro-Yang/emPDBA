def CG_percent(dna_fasta):
    num_CG = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'C' and seq[i + 1] == 'G':
                num_CG += 1
    per_CG = (num_CG/allpair) * 100
    return per_CG