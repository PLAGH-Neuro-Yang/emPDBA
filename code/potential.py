import struct
from math import sqrt
import os

aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'HIS': 'H', 'LYS': 'K', 'GLY': 'G',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W',
    'DA' :'DA', 'DT' :'DT', 'DG' :'DG', 'DC' :'DC'
}

dssp_codes = {
    'B': 'X', 'G': 'X', 'T': 'X',
    'H': 'Y', 'S': 'Y', ' ': 'Y',
    'E': 'Z', 'I': 'Z'
}

dna_location = {
    'DA': 2, 'DT': 3, 'DG': 4, 'DC': 5
}

def calcDist(p1, p2):
    tmp = pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2)
    tmp = sqrt(tmp)
    return tmp


def getInterfacepair(filename, chain1, chain2,dssp_data,pairfile):
    in_file = open(filename)
    pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
    A, B, results = [], [], []
    for line in in_file:
        if line[0:4] == "ATOM":
            col = struct.unpack(pdb_format, line.encode())
            a_name = col[3].strip().decode("utf-8")
            chain = col[7].strip().decode("utf-8")
            amino = col[5].strip().decode("utf-8")
            numer = col[8].strip().decode("utf-8")
            x = col[11].strip().decode("utf-8")
            y = col[12].strip().decode("utf-8")
            z = col[13].strip().decode("utf-8")
            if amino in aa_codes:
                amino = aa_codes[amino]
                if chain == chain1:
                    for i in range(len(dssp_data)):
                        if chain == dssp_data[i][2]:
                            if amino == dssp_data[i][0]:
                                if numer == dssp_data[i][1]:
                                    A.append([amino,numer,chain, x, y, z,dssp_data[i][3]])
                if chain == chain2:
                    if amino in aa_codes:
                        B.append([amino,numer,chain, x, y, z])

    for i in range(len(A)):
        for j in range(len(B)):
            v1 = (float(A[i][3]), float(A[i][4]), float(A[i][5]))
            v2 = (float(B[j][3]), float(B[j][4]), float(B[j][5]))
            tmp = calcDist(v1, v2)
            if tmp < 5:
                results.append((A[i][0],A[i][1],A[i][2],A[i][6],B[j][0],B[j][1],B[j][2]))

    result = list(set(results))
    result.sort(key=results.index)
    for line in result:
        for element in line:
            pairfile.write(str(element) + '\t')
        pairfile.write('\n')


def dssp(protein,dsspfile):
    dssp_path = './code/dssp/dsspcmbi'
    os.system('chmod 777 ' + dssp_path)
    os.system(dssp_path + ' ' + protein + ' > ' + dsspfile)

def calpotential(pairfile):
    deltaG = []
    deltaGall = 0
    for pair in pairfile:
        pair = pair.strip('\n').split('\t')
        amino = pair[0]
        dssp = dssp_codes[pair[3]]
        nucleic = pair[4]
        for aa_type in open("./code/60-4.data"):
            aa_type = aa_type.strip('\n').split('\t')
            amino_type = aa_type[0]
            dssp_type = aa_type[1]
            if amino == amino_type:
                if dssp == dssp_type:
                    deltaG.append(aa_type[dna_location[nucleic]])
    for g in range(0,len(deltaG)):
        deltaGall = deltaGall + float(deltaG[g])
    deltaGall = round(deltaGall,6)
    return deltaGall

def potential(complex,protein_chain,dna_chain,dsspfile):
    dssp(complex,dsspfile)
    dssp_data = []
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    line = open(dsspfile).readlines()
    title_index = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        num = col[1].strip().decode("utf-8")
        chain = col[2].strip().decode("utf-8")
        amino = col[3].strip().decode("utf-8")
        type = col[4].strip().decode("utf-8")
        if type == '':
            dssp_data.append([amino,num,chain, ' '])
        else:
            dssp_data.append([amino,num,chain, type])
    pairfile = open('./feature/pair.txt', 'w')
    for m in range(len(protein_chain)):
        for n in range(len(dna_chain)):
            getInterfacepair(complex, protein_chain[m],dna_chain[n],dssp_data,pairfile)
    pairfile.close()
    pairfile = open('./feature/pair.txt','r')
    pairwise = calpotential(pairfile)
    return pairwise

