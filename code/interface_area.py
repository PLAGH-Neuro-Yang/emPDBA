def interface_area(complex,protein,dna):
    import os
    tclname1 = './feature/sasa_complex.tcl'
    tclfile1 = open(tclname1,'w')
    tclfile1.write('mol load pdb ' + complex + '\n')
    tclfile1.write('set a [measure sasa 1.4 [atomselect top "all"]]\n')
    tclfile1.write('set outfile [open ./feature/sasa_complex.dat w]\n')
    tclfile1.write('puts $outfile $a\n')
    tclfile1.write('close $outfile\n')
    tclfile1.write('quit')
    tclfile1.close()
    os.system('vmd -dispdev text -e ./feature/sasa_complex.tcl')
    csasa = open('./feature/sasa_complex.dat').readlines()
    sasa_complex = float(csasa[0].split()[0])

    tclname2 = './feature/sasa_protein.tcl'
    tclfile2 = open(tclname2, 'w')
    tclfile2.write('mol load pdb ' + protein + '\n')
    tclfile2.write('set a [measure sasa 1.4 [atomselect top "all"]]\n')
    tclfile2.write('set outfile [open ./feature/sasa_protein.dat w]\n')
    tclfile2.write('puts $outfile $a\n')
    tclfile2.write('close $outfile\n')
    tclfile2.write('quit')
    tclfile2.close()
    os.system('vmd -dispdev text -e ./feature/sasa_protein.tcl')
    psasa = open('./feature/sasa_protein.dat').readlines()
    sasa_protein = float(psasa[0].split()[0])

    tclname3 = './feature/sasa_dna.tcl'
    tclfile3 = open(tclname3, 'w')
    tclfile3.write('mol load pdb ' + dna + '\n')
    tclfile3.write('set a [measure sasa 1.4 [atomselect top "all"]]\n')
    tclfile3.write('set outfile [open ./feature/sasa_dna.dat w]\n')
    tclfile3.write('puts $outfile $a\n')
    tclfile3.write('close $outfile\n')
    tclfile3.write('quit')
    tclfile3.close()
    os.system('vmd -dispdev text -e ./feature/sasa_dna.tcl')
    dsasa = open('./feature/sasa_dna.dat').readlines()
    sasa_dna = float(dsasa[0].split()[0])

    sasa = (sasa_protein+sasa_dna-sasa_complex)/2
    return sasa
