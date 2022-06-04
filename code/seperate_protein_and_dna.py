aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}

dna_codes = {
    'DA': 'A', 'DT': 'T', 'DG': 'G', 'DC': 'C',
}

def seperate_protein_and_dna(pdb,protein,dna):
    outpro = open(protein, 'w')
    outdna = open(dna, 'w')
    complex = open(pdb, 'r')
    for line in complex:
        line_term = line.split()
        if line_term[0] == 'ATOM':
            if line_term[3] in dna_codes:
                outdna.write(line)
            if line_term[3] in aa_codes:
                outpro.write(line)
    complex.close()
    outpro.close()
    outdna.close()
