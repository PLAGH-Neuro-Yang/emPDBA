def GC_num(dna_fasta):
    num_GC = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'G' and seq[i + 1] == 'C':
                num_GC += 1
    return num_GC