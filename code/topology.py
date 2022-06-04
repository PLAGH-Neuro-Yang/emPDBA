#!/usr/bin/env python
#    -*- coding:utf-8 –*-

"""
@author cc
@time 2017.12.12
"""
import os
import math
import numpy
import networkx as nx


# construct adjacent matrix
def create_matrix(number, amount):
    matrix = []
    for i in range(0, number):
        tmp = []
        for j in range(0, number):
            tmp.append(amount)
        matrix.append(tmp)
    return matrix

def pdbread(filename):
    pdb=open(filename, 'r')
    content=pdb.read()
    pdb.close
    line= content.strip().split("\n")
    List=[]
    node_list=[]
    aminoacid = ['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TRP','TYR','SER','THR','CYS','MET','ASN','GLN','ASP','GLU','HIS','LYS','ARG','PTR','MSE','SEP','TPO','AYA']
    nucleiotide = [' DA',' DT',' DG',' DC','5IU','BRU','8OG','5CM','PED','6MA','MRG','ATM','UMP','6MI','2DT','MA7','TED','DOC','N2G',
        ' TT','3DR','DDG','64T','5PY','6OG','EDC','DUZ','A2M','PE6','UPE','18Q','18M','1CC','OMC','OMG','2JU','5HC','G47','5FC','85Y',' GS','8DT','NR1','5NC','SOS','C2S']
    for i in range(0,len(line)):
        if line[i][0:4]=='ATOM' or line[i][0:6]=='HETATM':
            resname=line[i][17:20]
            if (resname in aminoacid and (line[i][13] == 'C' and line[i][14] == 'A')) or (resname in nucleiotide and line[i][13] == 'P'):
                x_coord = float(line[i][30:38])
                y_coord = float(line[i][38:46])
                z_coord = float(line[i][46:54])
                List.append([x_coord, y_coord, z_coord])
                node_list.append(resname)
                n=len(node_list)
    adjacent_matrix = create_matrix(n, 0)
    for i in range(len(List)):
        for j in range(len(List)):
            if i == j:
                continue
            else:
                dis = juli(List[i],List[j])
                if node_list[i] in aminoacid and node_list[j] in aminoacid:
                    cutoff = 7
                elif  node_list[i] in nucleiotide and node_list[j] in nucleiotide:
                    cutoff = 13
                elif (node_list[i] in aminoacid and node_list[j] in nucleiotide) or (node_list[i] in nucleiotide and node_list[j] in aminoacid):
                    cutoff=10
            if dis < cutoff:
                adjacent_matrix[i][j]=1
                adjacent_matrix[j][i]=1
    return  adjacent_matrix,node_list,n

def juli(x1,x2):
    juli = float(math.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2))
    return  juli
def get_filename(path,filetype):
    name =[]
    final_name = []
    for root,dirs,files in os.walk(path):
        for i in files:
            if filetype in i:
                name.append(i.replace(filetype,''))
    final_name = [item +'.pdb' for item in name]
    return final_name


def a_degrees_centrality(complex):
    aminoacid = ['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TRP','TYR','SER','THR','CYS','MET','ASN','GLN','ASP','GLU','HIS','LYS','ARG','PTR','MSE','SEP','TPO','AYA']
    nucleiotide = [' DA',' DT',' DG',' DC','5IU','BRU','8OG','5CM','PED','6MA','MRG','ATM','UMP','6MI','2DT','MA7','TED','DOC','N2G',
        ' TT','3DR','DDG','64T','5PY','6OG','EDC','DUZ','A2M','PE6','UPE','18Q','18M','1CC','OMC','OMG','2JU','5HC','G47','5FC','85Y',' GS','8DT','NR1','5NC','SOS','C2S']

    adjacent_matrix,node_list,n=pdbread(complex)
    G=nx.Graph()
    G=nx.from_numpy_matrix(numpy.matrix(adjacent_matrix))
    if (nx.is_connected(G)==True):
        print ('The graph is connected.')
        Graph_information=nx.info(G)
        print (Graph_information)
        degrees_centrality=nx.degree_centrality(G)
        average_degrees_centrality=0
        for key in range(0,len(degrees_centrality)):
            average_degrees_centrality=average_degrees_centrality+degrees_centrality[key]
        adc = average_degrees_centrality/n
        return adc
    else:
        print('The graph is unconnected.')

def average_closeness_centrality(complex):
    aminoacid = ['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TRP','TYR','SER','THR','CYS','MET','ASN','GLN','ASP','GLU','HIS','LYS','ARG','PTR','MSE','SEP','TPO','AYA']
    nucleiotide = [' DA',' DT',' DG',' DC','5IU','BRU','8OG','5CM','PED','6MA','MRG','ATM','UMP','6MI','2DT','MA7','TED','DOC','N2G',
        ' TT','3DR','DDG','64T','5PY','6OG','EDC','DUZ','A2M','PE6','UPE','18Q','18M','1CC','OMC','OMG','2JU','5HC','G47','5FC','85Y',' GS','8DT','NR1','5NC','SOS','C2S']

    adjacent_matrix,node_list,n=pdbread(complex)
    G=nx.Graph()
    G=nx.from_numpy_matrix(numpy.matrix(adjacent_matrix))
    if (nx.is_connected(G)==True):
        print ('The graph is connected.')
        Graph_information=nx.info(G)
        print (Graph_information)
        closenesses_centrality=nx.closeness_centrality(G)
        average_closeness_centrality=0
        for key in range(0,len(closenesses_centrality)):
            average_closeness_centrality=average_closeness_centrality+closenesses_centrality[key]
        acc = average_closeness_centrality/n
        return acc
    else:
        print('The graph is unconnected.')