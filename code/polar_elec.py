def protein_negative_percent(proteinlabel):
    residues=0
    all=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        all+=1
        if label[1] == "D" or label[1] == "E":
            residues += 1
    percent = (residues/all) * 100
    return percent

def protein_positive_num(proteinlabel):
    residues = 0
    labelfile = open(proteinlabel, 'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[1] == "K" or label[1] == "R" or label[1] == "H":
            residues += 1
    return residues

def protein_polar_uncharge_percent(proteinlabel):
    residues=0
    all=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        all+=1
        if label[1] == "T" or label[1] == "S" or label[1] == "C" \
                or label[1] == "N" or label[1] == "Q" or label[1] == "Y" \
                or label[1] == "G":
            residues += 1
    percent = (residues/all) * 100
    return percent

def protein_nonpolar_num(proteinlabel):
    residues = 0
    labelfile = open(proteinlabel, 'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[1] == "A" or label[1] == "F" or label[1] == "I" \
                or label[1] == "L" or label[1] == "M" or label[1] == "P" \
                or label[1] == "V" or label[1] == "W":
            residues += 1
    return residues