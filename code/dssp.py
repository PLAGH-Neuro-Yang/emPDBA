import os
import struct

def dssp(protein,dsspfile):
    dssp_path = './code/dssp/dsspcmbi'
    os.system('chmod 777 ' + dssp_path)
    os.system(dssp_path + ' ' + protein + ' > ' + dsspfile)

def dssp_B_percent(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile,'r')
    line = file.readlines()
    title_index = 0
    B_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    allresidue = len(line)-title_index-1
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'B':
            B_num += 1
    B_per = (B_num/allresidue) * 100
    return B_per

def dssp_I_num(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    I_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'I':
            I_num += 1
    return I_num

def dssp_I_percent(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    I_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    allresidue = len(line) - title_index - 1
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'I':
            I_num += 1
    I_per = (I_num / allresidue) * 100
    return I_per

def dssp_T_num(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    T_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'T':
            T_num += 1
    return T_num

def dssp_S_num(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    S_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'S':
            S_num += 1
    return S_num

def dssp_H_num(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    H_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'H':
            H_num += 1
    return H_num

def dssp_G_percent(dsspfile):
    dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
    file = open(dsspfile, 'r')
    line = file.readlines()
    title_index = 0
    G_num = 0
    for i in range(len(line)):
        if line[i][0:3].strip() == '#':
            title_index = i
    allresidue = len(line) - title_index - 1
    for i in range(title_index + 1, len(line)):
        col = struct.unpack(dssp_format, line[i].encode())
        type = col[4].strip().decode("utf-8")
        if type == 'G':
            G_num += 1
    G_per = (G_num / allresidue) * 100
    return G_per