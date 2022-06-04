import os

def vdw_elec(complex):
    name = complex.split('/')[-1]
    list = open('./code/elec-vdw/lis','w')
    list.write(name)
    list.close()
    os.system('chmod 777 ./code/elec-vdw/vdw.pl')
    os.system('./code/elec-vdw/vdw.pl ./code/elec-vdw/lis > ./code/elec-vdw/elec-vdw.dat')


def vdw_rep(complex):
    vdw_elec(complex)
    vdwr=0
    vdw = open('./code/elec-vdw/elec-vdw.dat', 'r')
    for line in vdw:
        line = line.strip('\n').split()
        vdwr = float(line[6])
    vdw.close()
    return vdwr

def elec_short_attr(complex):
    vdw_elec(complex)
    esa = 0
    vdw = open('./code/elec-vdw/elec-vdw.dat', 'r')
    for line in vdw:
        line = line.strip('\n').split()
        esa = float(line[2])
    vdw.close()
    return esa

def elec_short_rep(complex):
    vdw_elec(complex)
    esp = 0
    vdw = open('./code/elec-vdw/elec-vdw.dat', 'r')
    for line in vdw:
        line = line.strip('\n').split()
        esp = float(line[3])
    vdw.close()
    return esp