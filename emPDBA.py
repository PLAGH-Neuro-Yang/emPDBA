import os
import sys
import joblib
import struct
sys.path.append('./code/')
from seperate_protein_and_dna import seperate_protein_and_dna
from double_type import double_type
from vdw import vdw_attr
from hbonds import hbonds
from calculate_binding_sites import calc_sites
from get_sequence import get_fasta_from_label
from binding_sites_polar_uncharge import num_binding_sites_polar_uncharge
from binding_sites_polar_uncharge import per_binding_sites_polar_uncharge
from binding_sites_DA import per_binding_sites_DA
from dssp_G_percent import dssp_G_percent
from dssp_S_num import dssp_S_num
from beta_sheet import num_beta,weight_beta
from GA_num import GA_num
from CT_num import CT_num
from CG_num import CG_num
from GG_num import GG_num
from AC_num import AC_num
from TA_num import TA_num
from GC_num import GC_num
from CG_percent import CG_percent
sys.path.append('./code/model/')
from ADA_model import ADA_model
from BAG_model import BAG_model
from DTR_model import DTR_model
from GBR_model import GBR_model
from RFR_model import RFR_model
from XGB_model import XGB_model

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
    x1 = (GA_num(dna_fasta)/5)*2-1
    x2 = (CG_percent(dna_fasta)/40)*2-1
    x3 = (per_binding_sites_polar_uncharge(protein_label)/63.15789474)*2-1
    x = [x1,x2,x3]
    return x

def feature_extraction_DII(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    dsspfile = './feature/protein.dssp'
    x1 = (dssp_G_percent(protein,dsspfile)/9.763313609)*2-1
    x2 = (hbonds(pdb)/74)*2-1
    x3 = (CT_num(dna_fasta)/7)*2-1
    x4 = (CG_num(dna_fasta)/6)*2-1
    x = [x1,x2,x3,x4]
    return x

def feature_extraction_DIII(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    x1 = (GG_num(dna_fasta)/21)*2-1
    x2 = ((per_binding_sites_DA(dna_label)-6.25)/(54.54545455-6.25))*2-1
    x3 = (vdw_attr(complex)/13.2491554)*2-1
    x = [x1,x2,x3]
    return x

def feature_extraction_M(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    dsspfile = './feature/protein.dssp'
    x1 = ((dssp_S_num(protein,dsspfile)-1)/106)*2-1
    x2 = (AC_num(dna_fasta)/4)*2-1
    x3 = (TA_num(dna_fasta)/2)*2-1
    x4 = (GG_num(dna_fasta)/6)*2-1
    x5 = (GC_num(dna_fasta)/8)*2-1
    x6 = (num_binding_sites_polar_uncharge(protein_label)/20)*2-1
    x = [x1,x2,x3,x4,x5,x6]
    return x

def feature_extraction_S(complex,protein,dna):
    protein_label = './feature/protein_label.data'
    dna_label = './feature/dna_label.data'
    calc_sites(protein, dna, protein_label, dna_label)
    dna_fasta = './feature/dna.fasta'
    get_fasta_from_label(dna_label, dna_fasta)
    beta_num,beta_weight = 0,0
    for chain in protein_chain:
        betafile = './feature/label_beta_' + chain + '.data'
        beta_num += num_beta(protein_label,protein,betafile)
        beta_weight += weight_beta(protein_label,protein,betafile)
    x1 = GC_num(dna_fasta)*2-1
    x2 = ((beta_num-6)/15)*2-1
    x3 = ((beta_weight-5819.33)/(15911.99-5819.33))*2-1
    x = [x1,x2,x3]
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
    elif type == 'S':
        x = feature_extraction_S(complex,protein,dna)
    else:
        print('Please input a correct complex type.')
    print('***** Feature extraction completed! *****')
    # print(x)

    #prediction
    print('***** Start prediction ... *****')
    y_predict = 0
    if type == 'DI':
        train = './code/model/DI.data'
        params_dtr = {'random_state': 67}
        params_xgb = {'n_estimators': 20, 'max_depth': 8, 'learning_rate': 0.28}
        params_rfr = {'n_estimators': 1, 'random_state': 21, 'criterion': 'mse'}
        params_ada = {'n_estimators': 228, 'random_state': 46, 'learning_rate': 0.9}
        params_bag = {'n_estimators': 2, 'random_state': 59}
        params_gbr = {'n_estimators': 27, 'max_depth': 3, 'min_samples_split': 2,'learning_rate': 0.4, 'loss': 'huber',
                      'random_state': 29}
        y_ensemble = [[DTR_model(train,x,params_dtr),XGB_model(train,x,params_xgb),
                      RFR_model(train,x,params_rfr),ADA_model(train,x,params_ada),
                      BAG_model(train,x,params_bag),GBR_model(train,x,params_gbr)]]
        model = joblib.load('./code/model/DI.m')
        y_predict = model.predict(y_ensemble)
    elif type == 'DII':
        train = './code/model/DII.data'
        params_dtr = {'random_state': 97}
        params_xgb = {'n_estimators': 88, 'max_depth': 5, 'learning_rate': 0.26}
        params_rfr = {'n_estimators': 14, 'random_state': 34, 'criterion': 'mse'}
        params_ada = {'n_estimators': 72, 'random_state': 79, 'learning_rate': 0.6}
        params_bag = {'n_estimators': 53, 'random_state': 34}
        params_gbr = {'n_estimators': 23, 'max_depth': 3, 'min_samples_split': 2,
                              'learning_rate': 0.4, 'loss': 'huber', 'random_state': 81}
        y_ensemble = [[DTR_model(train, x, params_dtr), XGB_model(train, x, params_xgb),
                      RFR_model(train, x, params_rfr), ADA_model(train, x, params_ada),
                      BAG_model(train, x, params_bag), GBR_model(train, x, params_gbr)]]
        model = joblib.load('./code/model/DII.m')
        y_predict = model.predict(y_ensemble)
    elif type == 'DIII':
        train = './code/model/DIII.data'
        params_dtr = {'random_state': 5}
        params_xgb = {'n_estimators': 18, 'max_depth': 6,'learning_rate': 0.9}
        params_rfr = {'n_estimators': 1, 'random_state': 97, 'criterion': 'mse'}
        params_ada = {'n_estimators': 16, 'random_state': 34,'learning_rate':0.97}
        params_bag = {'n_estimators': 1, 'random_state': 97}
        params_gbr = {'n_estimators': 3, 'max_depth': 3, 'min_samples_split': 2,
                              'learning_rate': 0.69, 'loss': 'huber', 'random_state': 4}
        y_ensemble = [[DTR_model(train, x, params_dtr), XGB_model(train, x, params_xgb),
                      RFR_model(train, x, params_rfr), ADA_model(train, x, params_ada),
                      BAG_model(train, x, params_bag), GBR_model(train, x, params_gbr)]]
        model = joblib.load('./code/model/DIII.m')
        y_predict = model.predict(y_ensemble)
    elif type == 'M':
        train = './code/model/M.data'
        params_dtr = {'random_state': 5}
        params_xgb = {'n_estimators': 45, 'max_depth': 3, 'learning_rate': 0.33}
        params_rfr = {'n_estimators': 19, 'random_state': 39, 'criterion': 'mse'}
        params_ada = {'n_estimators': 20, 'random_state': 71,'learning_rate':0.18}
        params_bag = {'n_estimators': 2, 'random_state': 27}
        params_gbr = {'n_estimators': 5, 'max_depth': 3, 'min_samples_split': 2,
                              'learning_rate': 0.5, 'loss': 'huber', 'random_state': 78}
        y_ensemble = [[DTR_model(train, x, params_dtr), XGB_model(train, x, params_xgb),
                      RFR_model(train, x, params_rfr), ADA_model(train, x, params_ada),
                      BAG_model(train, x, params_bag), GBR_model(train, x, params_gbr)]]
        model = joblib.load('./code/model/M.m')
        y_predict = model.predict(y_ensemble)
    elif type == 'S':
        train = './code/model/S.data'
        params_dtr = {'random_state': 11}
        params_xgb = {'n_estimators': 30, 'max_depth': 5,'learning_rate': 0.32}
        params_rfr = {'n_estimators': 1, 'random_state': 44, 'criterion': 'mse'}
        params_ada = {'n_estimators': 1, 'random_state': 73,'learning_rate':0.1}
        params_bag = {'n_estimators': 1, 'random_state': 9}
        params_gbr = {'n_estimators': 21, 'max_depth': 3, 'min_samples_split': 2,
                              'learning_rate': 0.9, 'loss': 'huber', 'random_state': 77}
        y_ensemble = [[DTR_model(train, x, params_dtr), XGB_model(train, x, params_xgb),
                      RFR_model(train, x, params_rfr), ADA_model(train, x, params_ada),
                      BAG_model(train, x, params_bag), GBR_model(train, x, params_gbr)]]
        model = joblib.load('./code/model/S.m')
        y_predict = model.predict(y_ensemble)
    print('***** Prediction completed! *****')
    print('The predicted binding affinity of ' + pdbname + ' is -' + str("%.2f" % y_predict) + 'kcal/mol.')



