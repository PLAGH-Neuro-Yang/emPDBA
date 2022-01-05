def TA_num(dna_fasta):
    num_TA = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'T' and seq[i + 1] == 'A':
                num_TA += 1
    return num_TA