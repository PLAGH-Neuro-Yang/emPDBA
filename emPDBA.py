import os
import sys
import joblib
import struct
import numpy as np
sys.path.append('./code/')
from seperate_protein_and_dna import seperate_protein_and_dna
from double_type import double_type
from calculate_binding_sites import calc_sites
from get_sequence import get_fasta_from_label
from polar_elec import protein_negative_percent,protein_positive_num,protein_polar_uncharge_percent,protein_nonpolar_num
from dssp import dssp,dssp_B_percent,dssp_I_num,dssp_I_percent,dssp_T_num,dssp_S_num,dssp_H_num,dssp_G_percent
from binding_sites import per_binding_sites_positive,per_binding_sites_polar_uncharge,num_binding_sites_DT,\
    per_binding_sites_nonpolar,per_binding_sites_negative,per_binding_sites_DA,per_binding_sites_DC
from interface_area import interface_area
from vdw_elec import vdw_rep,elec_short_attr,elec_short_rep
from topology import a_degrees_centrality,average_closeness_centrality
from volume import volumeDNA,volumecomplex
from reaction_field_energy import reaction_field_energy_DNA,change_reaction_field_energy
from alpha_helix import num_alpha,weight_alpha
from mass import massDNA
from hbonds import hbonds
from potential import potential
from dinucleotide import CA_num,CG_num,AC_num,AA_num,AT_num,AT_percent,TT_num,TC_num,\
    TA_percent,TA_num,GA_num,GA_percent,TG_num,CA_percent,GC_percent,GC_num


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

def feature_extraction_DI(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    dsspfile = './feature/protein.dssp'
    dssp(protein,dsspfile)
    x1 = ((protein_negative_percent(protein_label)-9.045)/(17.318-9.045))*2-1
    x2 = (dssp_B_percent(dsspfile)/3.135889)*2-1
    x3 = (dssp_I_num(dsspfile)/6)*2-1
    x4 = (dssp_I_percent(dsspfile)/0.421644)*2-1
    x5 = ((dssp_T_num(dsspfile)-14)/(200-14))*2-1
    x6 = ((dssp_S_num(dsspfile)-9)/(242-9))*2-1
    x7 = (TA_num(dna_fasta)/18)*2-1
    x8 = (GC_num(dna_fasta)/9)*2-1
    x9 = (CA_num(dna_fasta)/6)*2-1
    x10 = ((per_binding_sites_positive(protein_label)-9.090909)/(83.33333-9.090909))*2-1
    x11 = ((interface_area(complex,protein,dna)-184.3369)/(2690.508-184.3369))*2-1
    x12 = (vdw_rep(complex)/0.332572)*2-1
    x13 = ((a_degrees_centrality(pdb)-0.00431)/(0.034707-0.00431))*2-1
    x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13]
    return x

def feature_extraction_DII(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    volumefile = './feature/volume.data'
    reaction_field_energy_file = './feature/reaction_field_energy.data'
    x1 = ((protein_positive_num(protein_label)-13)/(142-13))*2-1
    x2 = (AA_num(dna_fasta)/14)*2-1
    x3 = (AT_num(dna_fasta)/18)*2-1
    x4 = (AT_percent(dna_fasta)/30)*2-1
    x5 = (TA_num(dna_fasta)/16)*2-1
    x6 = (TT_num(dna_fasta)/14)*2-1
    x7 = (CG_num(dna_fasta)/12)*2-1
    x8 = ((per_binding_sites_polar_uncharge(protein_label)-13.63636)/(69.1358-13.63636))*2-1
    x9 = (num_binding_sites_DT(dna_label)/37)*2-1
    x10 = ((interface_area(complex,protein,dna)-503.8757)/(6960.962-503.8757))*2-1
    x11 = ((elec_short_attr(complex)+444.096)/(-29.7624+444.096))*2-1
    x12 = (vdw_rep(complex)/61.49715)*2-1
    x13 = ((volumeDNA(volumefile)-3615.49)/(34412.5-3615.49))*2-1
    x14 = ((reaction_field_energy_DNA(reaction_field_energy_file)-11962.69)/(357171.6-11962.69))*2-1
    x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14]
    return x

def feature_extraction_DIII(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    dsspfile = './feature/protein.dssp'
    dssp(protein, dsspfile)
    reaction_field_energy_file = './feature/reaction_field_energy.data'
    alpha_num, alpha_weight = 0, 0
    for chain in protein_chain:
        alphafile = './feature/label_alpha_' + chain + '.data'
        alpha_num += num_alpha(protein_label, protein, alphafile)
        alpha_weight += weight_alpha(protein_label, protein, alphafile)
    x1 = ((protein_polar_uncharge_percent(protein_label)-17.46)/(57.692-17.46))*2-1
    x2 = (dssp_H_num(dsspfile) / 505) * 2 - 1
    x3 = (alpha_num / 35) * 2 - 1
    x4 = (alpha_weight / 66797.21) * 2 - 1
    x5 = (AT_percent(dna_fasta) / 33.33333) * 2 - 1
    x6 = (TA_percent(dna_fasta) / 33.33333) * 2 - 1
    x7 = (TC_num(dna_fasta) / 16) * 2 - 1
    x8 = (GA_num(dna_fasta) / 16) * 2 - 1
    x9 = (GC_percent(dna_fasta) / 30) * 2 - 1
    x10 = (CA_percent(dna_fasta) / 18.18182) * 2 - 1
    x11 = ((massDNA(dna) - 3723.4) / (87240.9 - 3723.4)) * 2 - 1
    x12 = (hbonds(pdb) / 154) * 2 - 1
    x13 = ((potential(pdb, protein_chain, dna_chain, dsspfile) + 106.862) / (10.655 + 106.862)) * 2 - 1
    x14 = (vdw_rep(complex) / 10.23937) * 2 - 1
    x15 = ((reaction_field_energy_DNA(reaction_field_energy_file) - 19028.04) / (2313147 - 19028.04)) * 2 - 1
    x16 = ((change_reaction_field_energy(reaction_field_energy_file) + 334075) / (2702807 + 334075)) * 2 - 1
    x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16]
    return x

def feature_extraction_M(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    dsspfile = './feature/protein.dssp'
    dssp(protein, dsspfile)
    volumefile = './feature/volume.data'
    x1 = ((protein_nonpolar_num(protein_label) - 8) / (706 - 8)) * 2 - 1
    x2 = (dssp_G_percent(dsspfile) / 11.146) * 2 - 1
    x3 = (AC_num(dna_fasta) / 8) * 2 - 1
    x4 = (TA_percent(dna_fasta) / 33.33333) * 2 - 1
    x5 = (TG_num(dna_fasta) / 5) * 2 - 1
    x6 = (GA_percent(dna_fasta) / 28.57143) * 2 - 1
    x7 = (CA_percent(dna_fasta) / 25) * 2 - 1
    x8 = (per_binding_sites_nonpolar(protein_label) / 100) * 2 - 1
    x9 = (per_binding_sites_negative(protein_label) / 19) * 2 - 1
    x10 = (per_binding_sites_DA(dna_label) / 100) * 2 - 1
    x11 = (per_binding_sites_DC(dna_label) / 100) * 2 - 1
    x12 = ((elec_short_rep(complex) - 1.154123) / (1322.437 - 1.154123)) * 2 - 1
    x13 = (vdw_rep(complex) / 81.32054) * 2 - 1
    x14 = ((average_closeness_centrality(complex) - 0.079952) / (0.353193 - 0.079952)) * 2 - 1
    x15 = ((volumecomplex(volumefile) - 7429.38) / (274073 - 7429.38)) * 2 - 1
    x = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15]
    return x


if __name__ == '__main__':
    pdb = os.sys.argv[1]
    type = os.sys.argv[2]
    pdbname = pdb.split('/')[-1].split('.')[0]

    #get complex.pdb, protein.pdb, dna.pdb and chain_id
    makecomplexfile = open('./feature/complex.pdb','w')
    protein_chain = ''
    dna_chain = ''
    pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s6s4s2s3s'
    for line in open(pdb):
        if line[0:4] == "ATOM":
            col = struct.unpack(pdb_format, line.encode())
            if col[18].strip().decode("utf-8") != 'H':
                makecomplexfile.write(line)
            res_name = col[5].strip().decode("utf-8")
            chain = col[7].strip().decode("utf-8")
            if res_name in aa_codes:
                if chain not in protein_chain:
                    protein_chain += chain
            if res_name in dna_codes:
                if chain not in dna_chain:
                    dna_chain += chain
    makecomplexfile.close()
    complex = './feature/complex.pdb'
    protein = './feature/protein.pdb'
    dna = './feature/dna.pdb'
    seperate_protein_and_dna(complex,protein,dna)
    if type == 'D':
        type = double_type(protein)

    #feature extraction
    print('***** Start feature extraction ... *****')
    x = []
    if type == 'DI':
        x = feature_extraction_DI(complex,protein,dna)
    elif type == 'DII':
        x = feature_extraction_DII(complex,protein,dna)
    elif type == 'DIII':
        x = feature_extraction_DIII(complex,protein,dna)
    elif type == 'M':
        x = feature_extraction_M(complex,protein,dna)
    else:
        print('Please input a correct complex type.')
    print('***** Feature extraction completed! *****')
    # print(x)

    #prediction
    print('***** Start prediction ... *****')
    y_predict = 0
    if type == 'DI':
        model = joblib.load('./code/model/DI.m')
        y_predict = model.predict(np.array([x]))
    elif type == 'DII':
        model = joblib.load('./code/model/DII.m')
        y_predict = model.predict(np.array([x]))
    elif type == 'DIII':
        model = joblib.load('./code/model/DIII.m')
        y_predict = model.predict(np.array([x]))
    elif type == 'M':
        model = joblib.load('./code/model/M.m')
        y_predict = model.predict(np.array([x]))
    print('***** Prediction completed! *****')
    print('The predicted binding affinity of ' + pdbname + ' is ' + str("%.2f" % y_predict) + 'kcal/mol.')



