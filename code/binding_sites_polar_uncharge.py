def num_binding_sites_polar_uncharge(proteinlabel):
    residues=0
    labelfile = open(proteinlabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            if label[1] == "T" or label[1] == "S" or label[1] == "C" \
                    or label[1] == "N" or label[1] == "Q" or label[1] == "Y" \
                    or label[1] == "G":
                residues += 1
    return residues

def per_binding_sites_polar_uncharge(proteinlabel):
    num = num_binding_sites_polar_uncharge(proteinlabel)
    sites = 0
    labelfile = open(proteinlabel, 'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            sites += 1
    per = (num/sites) * 100
    return per