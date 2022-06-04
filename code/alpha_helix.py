import struct
import os

aa_codes = {
    'A': 89.09,  'C': 121.16, 'D': 133.10, 'E': 147.13,
    'F': 165.19, 'G': 75.07,  'H': 155.16, 'K': 146.19,
    'I': 131.17, 'L': 131.17, 'M': 149.21, 'N': 132.12,
    'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
    'T': 119.16, 'V': 117.15, 'Y': 181.19, 'W': 204.22
}

def label_alpha(labelfile, protein,outfile):
    dssp_path = './code/dssp/dsspcmbi'
    dsspfile = './feature/protein.dssp'
    os.system('chmod 777 ' + dssp_path)
    os.system(dssp_path + ' ' + protein + ' > ' + dsspfile)
    labelfile = open(labelfile, 'r')
    label = labelfile.readlines()
    pre_chain = []
    dsspinput = open(dsspfile + "", 'r')
    dssp = dsspinput.readlines()
    for eachlabel in label:
        if not len(eachlabel.strip()):
            continue
        onelabel = eachlabel.strip('\n').split()
        chain = onelabel[0]
        amino = onelabel[1] + onelabel[2]
        if chain != pre_chain:
            outfile = open('./feature/label_alpha_' + chain + '.data', 'w')
        dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
        title_index = 0
        for i in range(len(dssp)):
            if dssp[i][0:3].strip() == '#':
                title_index = i  # 确定dssp标题的行index
        for i in range(title_index + 1, len(dssp)):
            col = struct.unpack(dssp_format, dssp[i].encode())
            chainid = col[2].strip().decode("utf-8")
            type = col[4].strip().decode("utf-8")
            aminonum = col[1].strip().decode("utf-8")
            aminoname = col[3].strip().decode("utf-8")
            aminoid = aminoname+aminonum
            if chainid == chain:
                if aminoid == amino:
                    if 'H' in type:
                        outfile.write(aminoname + '\t' + aminonum + '\t1\n')
                    else:
                        outfile.write(aminoname + '\t' + aminonum + '\t0\n')
        pre_chain = chain
    labelfile.close()
    dsspinput.close()
    outfile.close()

def num_alpha(labelfile,dsspfile,alphafile):
    label_alpha(labelfile, dsspfile,alphafile)
    alphalabel = open(alphafile, 'r')
    label = alphalabel.readlines()
    alpha=[]
    num=0
    for i in range(len(label)-1):
        alphai = label[i].strip('\n').split('\t')[2]
        alphai1 = label[i+1].strip('\n').split('\t')[2]
        if alphai=='0' and alphai1=='0':
            continue
        elif alphai=='0' and alphai1=='1':
            num=num+1
        elif alphai=='1' and alphai1=='1':
            num=num+1
        elif alphai=='1' and alphai1=='0':
            alpha.append(num)
            num=0
        if i==len(label)-2 and alphai1=='1':
            alpha.append(num)
            num = 0
    alpha_num = len(alpha)
    alphalabel.close()
    return alpha_num

def weight_alpha(labelfile,dsspfile,alphafile):
    label_alpha(labelfile, dsspfile, alphafile)
    alphalabel = open(alphafile, 'r')
    label = alphalabel.readlines()
    alphaweight = []
    for i in range(len(label)):
        alphai = label[i].strip('\n').split('\t')[2]
        amino = label[i].strip('\n').split('\t')[0]
        if alphai=='1':
            alphaweight.append(aa_codes[amino])
    alpha_weight = sum(alphaweight)
    alphalabel.close()
    return alpha_weight