def GA_num(dna_fasta):
    num_GA = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'G' and seq[i + 1] == 'A':
                num_GA += 1
    return num_GA