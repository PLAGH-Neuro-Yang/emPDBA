from math import sqrt
import struct

def calc_dist(p1, p2):
    tmp = pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2)
    tmp = sqrt(tmp)
    return tmp

def double_type(protein):
    dnafile = open('./feature/dna.pdb', 'r')
    dnainput = dnafile.readlines()
    output = []
    with open(protein, 'r') as proteininput:
        resname = []
        resnum = []
        label = "unbind"
        pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
        for proteinline in proteininput:
            col_pro = struct.unpack(pdb_format, proteinline.encode())
            aaname = col_pro[5].strip().decode("utf-8")
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
                output.append(out_label)
                label = 'unbind'
            resname = aaname
            resnum = aanum
            for dnaline in dnainput:
                col_dna = struct.unpack(pdb_format, dnaline.encode())
                dnax = col_dna[11].strip().decode("utf-8")
                dnay = col_dna[12].strip().decode("utf-8")
                dnaz = col_dna[13].strip().decode("utf-8")
                aacor = [float(aax), float(aay), float(aaz)]
                dnacor = [float(dnax), float(dnay), float(dnaz)]
                dist = calc_dist(aacor, dnacor)
                if dist - 5 <= 0.0:
                    label = "bind"
                    break
        out_label = '-'
        if label == 'bind':
            out_label = '+'
        output.append(out_label)
    dnafile.close()
    rate = (output.count('+') / len(output)) * 100
    if rate <= 10:
        return 'DI'
    elif rate >= 20:
        return 'DIII'
    else:
        return 'DII'