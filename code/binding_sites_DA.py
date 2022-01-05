def num_binding_sites_DA(dnalabel):
    residues=0
    labelfile = open(dnalabel,'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            if label[1] == "A" :
                residues += 1
    return residues

def per_binding_sites_DA(dnalabel):
    num = num_binding_sites_DA(dnalabel)
    sites = 0
    labelfile = open(dnalabel, 'r')
    for label in labelfile:
        label = label.strip('\n').split()
        if label[3] == '+':
            sites += 1
    per = (num/sites) * 100
    return per