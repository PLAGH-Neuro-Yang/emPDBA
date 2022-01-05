def hbonds(pdb):
    import os
    tclname = './feature/hbonds.tcl'
    tclfile = open(tclname,'w')
    tclfile.write('mol load pdb ' + pdb + '\n')
    tclfile.write('set sel1 [atomselect 0 protein]\n')
    tclfile.write('set sel2 [atomselect 0 nucleic]\n')
    tclfile.write('set a [llength [lindex [measure hbonds 3.5 45 $sel1 $sel2] 0]]\n')
    tclfile.write('set b [llength [lindex [measure hbonds 3.5 45 $sel2 $sel1] 0]]\n')
    tclfile.write('set c [expr $a+$b]\n')
    tclfile.write('set outfile [open ./feature/hbonds.dat w]\n')
    tclfile.write('puts $outfile $c\n')
    tclfile.write('close $outfile\n')
    tclfile.write('quit')
    tclfile.close()
    os.system('vmd -dispdev text -e ./feature/hbonds.tcl')
    hbondsnum = open('./feature/hbonds.dat').readlines()
    num = float(hbondsnum[0].split()[0])
    return num
