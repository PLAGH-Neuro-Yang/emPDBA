def CT_num(dna_fasta):
    num_CT = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'C' and seq[i + 1] == 'T':
                num_CT += 1
    return num_CT