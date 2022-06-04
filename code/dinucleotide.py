def CA_num(dna_fasta):
    num_CA = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'C' and seq[i + 1] == 'A':
                num_CA += 1
    return num_CA

def CG_num(dna_fasta):
    num_CG = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'C' and seq[i + 1] == 'G':
                num_CG += 1
    return num_CG

def AC_num(dna_fasta):
    num_AC = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'A' and seq[i + 1] == 'C':
                num_AC += 1
    return num_AC

def AA_num(dna_fasta):
    num_AA = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'A' and seq[i + 1] == 'A':
                num_AA += 1
    return num_AA

def AT_num(dna_fasta):
    num_AT = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'A' and seq[i + 1] == 'T':
                num_AT += 1
    return num_AT

def AT_percent(dna_fasta):
    num_AT = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'A' and seq[i + 1] == 'T':
                num_AT += 1
    per_AT = (num_AT/allpair) * 100
    return per_AT

def TT_num(dna_fasta):
    num_TT = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'T' and seq[i + 1] == 'T':
                num_TT += 1
    return num_TT

def TC_num(dna_fasta):
    num_TC = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'T' and seq[i + 1] == 'C':
                num_TC += 1
    return num_TC

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

def TA_percent(dna_fasta):
    num_TA = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'T' and seq[i + 1] == 'A':
                num_TA += 1
    per_TA = (num_TA/allpair) * 100
    return per_TA

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

def GA_percent(dna_fasta):
    num_GA = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'G' and seq[i + 1] == 'A':
                num_GA += 1
    per_GA = (num_GA/allpair) * 100
    return per_GA

def TG_num(dna_fasta):
    num_TG = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        for i in range(len(seq) - 1):
            if seq[i] == 'T' and seq[i + 1] == 'G':
                num_TG += 1
    return num_TG

def CA_percent(dna_fasta):
    num_CA = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'C' and seq[i + 1] == 'A':
                num_CA += 1
    per_CA = (num_CA/allpair) * 100
    return per_CA

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

def GC_percent(dna_fasta):
    num_GC = 0
    allpair = 0
    fasta = open(dna_fasta, 'r')
    for eachline in fasta:
        if '>' in eachline:
            continue
        seq = eachline.strip()
        allpair += len(seq)-1
        for i in range(len(seq) - 1):
            if seq[i] == 'G' and seq[i + 1] == 'C':
                num_GC += 1
    per_GC = (num_GC/allpair) * 100
    return per_GC