from math import sqrt
import struct

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

def calc_dist(p1, p2):
    tmp = pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2)
    tmp = sqrt(tmp)
    return tmp

# calculate and label the binding sites in protein
def label_protein_binding_sites(protein,dna):
    dnafile = open(dna, 'r')
    dnainput = dnafile.readlines()
    output = ''
    with open(protein, 'r') as proteininput:
        resname = ''
        resnum = ''
        chain0 = ''
        label = "unbind"
        pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
        for proteinline in proteininput:
            col_pro = struct.unpack(pdb_format, proteinline.encode())
            aaname = col_pro[5].strip().decode("utf-8")
            chain1 = col_pro[7].strip().decode("utf-8")
            aanum = col_pro[8].strip().decode("utf-8")
            aax = col_pro[11].strip().decode("utf-8")
            aay = col_pro[12].strip().decode("utf-8")
            aaz = col_pro[13].strip().decode("utf-8")
            if label == 'bind' and aaname == resname and aanum == resnum:
                continue
            if (aaname != resname or aanum != resnum) and len(resname):
                out_label = '-'
                if label == 'bind':
                    out_label = '+'
                reslabel = chain0 + '\t' + aa_codes[resname] + '\t' + resnum + '\t' + out_label + '\n'
                output += reslabel
                label = 'unbind'
            resname = aaname
            resnum = aanum
            chain0 = chain1
            for dnaline in dnainput:
                col_dna = struct.unpack(pdb_format, dnaline.encode())
                dnax = col_dna[11].strip().decode("utf-8")
                dnay = col_dna[12].strip().decode("utf-8")
                dnaz = col_dna[13].strip().decode("utf-8")
                aacor = [float(aax),float(aay),float(aaz)]
                dnacor = [float(dnax), float(dnay), float(dnaz)]
                dist = calc_dist(aacor,dnacor)
                if dist - 4 <= 0.0:
                    label = "bind"
                    break
        out_label = '-'
        if label == 'bind':
            out_label = '+'
        reslabel = chain0 + '\t' + aa_codes[resname] + '\t' + resnum + '\t' + out_label + '\n'
        output += reslabel
    dnafile.close()
    return output


# calculate and label the binding sites in dna
def label_dna_binding_sites(protein,dna):
    proteinfile = open(protein, 'r')
    proteininput = proteinfile.readlines()
    output = ''
    with open(dna, 'r') as dnainput:
        resname = ''
        resnum = ''
        chain0 = ''
        label = "unbind"
        pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
        for dnaline in dnainput:
            col_dna = struct.unpack(pdb_format, dnaline.encode())
            dnaname = col_dna[5].strip().decode("utf-8")
            chain1 = col_dna[7].strip().decode("utf-8")
            dnanum = col_dna[8].strip().decode("utf-8")
            dnax = col_dna[11].strip().decode("utf-8")
            dnay = col_dna[12].strip().decode("utf-8")
            dnaz = col_dna[13].strip().decode("utf-8")
            if label == 'bind' and dnaname == resname and dnanum == resnum:
                continue
            if (dnaname != resname or dnanum != resnum) and len(resname):
                out_label = '-'
                if label == 'bind':
                    out_label = '+'
                reslabel = chain0 + '\t' + dna_codes[resname] + '\t' + resnum + '\t' + out_label + '\n'
                output += reslabel
                label = 'unbind'
            resname = dnaname
            resnum = dnanum
            chain0 = chain1
            for proteinline in proteininput:
                col_pro = struct.unpack(pdb_format, proteinline.encode())
                aax = col_pro[11].strip().decode("utf-8")
                aay = col_pro[12].strip().decode("utf-8")
                aaz = col_pro[13].strip().decode("utf-8")
                aacor = [float(aax), float(aay), float(aaz)]
                dnacor = [float(dnax), float(dnay), float(dnaz)]
                dist = calc_dist(aacor,dnacor)
                if dist - 4 <= 0.0:
                    label = "bind"
                    break
        out_label = '-'
        if label == 'bind':
            out_label = '+'
        reslabel = chain0 + '\t' + dna_codes[resname] + '\t' + resnum + '\t' + out_label + '\n'
        output += reslabel
    proteinfile.close()
    return output


# output the label of binding site residues
def calc_sites(protein,dna,protein_label,dna_label):
    proteinlabel = label_protein_binding_sites(protein,dna)
    dnalabel = label_dna_binding_sites(protein,dna)
    outfile_protein = open(protein_label,'w')
    outfile_protein.write(proteinlabel)
    outfile_dna = open(dna_label,'w')
    outfile_dna.write(dnalabel)
    outfile_protein.close()
    outfile_dna.close()


