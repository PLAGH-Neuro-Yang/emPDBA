import os

def vdw_attr(complex):
    name = complex.split('/')[-1]
    list = open('./code/elec-vdw/lis','w')
    list.write(name)
    list.close()
    os.system('chmod 777 ./code/elec-vdw/vdw.pl')
    os.system('./code/elec-vdw/vdw.pl ./code/elec-vdw/lis > ./code/elec-vdw/elec-vdw.dat')
    vdwa = 0
    vdw = open('./code/elec-vdw/elec-vdw.dat','r')
    for line in vdw:
        line = line.strip('\n').split()
        vdwa = float(line[6])
    vdw.close()
    os.system('rm ./code/elec-vdw/elec-vdw.dat')
    return vdwa