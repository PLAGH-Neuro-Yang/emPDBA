def per_binding_sites_positive(proteinlabel):
    residues=0
    sites=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            sites+=1
            if label[1] == "K" or label[1] == "R" or label[1] == "H":
                residues += 1
    per = (residues / sites) * 100
    return per

def per_binding_sites_polar_uncharge(proteinlabel):
    residues=0
    all=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            all+=1
            if label[1] == "T" or label[1] == "S" or label[1] == "C" \
                    or label[1] == "N" or label[1] == "Q" or label[1] == "Y" \
                    or label[1] == "G":
                residues += 1
    per = (residues/all)*100
    return per

def num_binding_sites_DT(dnalabel):
    residues=0
    labelfile = open(dnalabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            if label[1] == "T":
                residues += 1
    return residues

def per_binding_sites_nonpolar(proteinlabel):
    residues=0
    all=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            all+=1
            if label[1] == "A" or label[1] == "F" or label[1] == "I" \
                    or label[1] == "L" or label[1] == "M" or label[1] == "P" \
                    or label[1] == "V" or label[1] == "W":
                residues += 1
    per = (residues/all)*100
    return per

def per_binding_sites_negative(proteinlabel):
    residues=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            if label[1] == "D" or label[1] == "E":
                residues += 1
    return residues

def per_binding_sites_DA(dnalabel):
    residues=0
    all =0
    labelfile = open(dnalabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            all+=1
            if label[1] == "A":
                residues += 1
    per = (residues/all)*100
    return per

def per_binding_sites_DC(dnalabel):
    residues=0
    all =0
    labelfile = open(dnalabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            all+=1
            if label[1] == "C":
                residues += 1
    per = (residues/all)*100
    return per