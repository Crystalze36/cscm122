#!/usr/bin/env python
# coding: utf-8

# In[2]:


from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None

def preProcess(datas, k, thrd):
    dkm = defaultdict(int)
    pkm = len(datas[0][0]) - k
    for data in datas:
        i = 0
        while i < pkm + 1:
            dkm[data[0][i : i + k + 1]], dkm[data[1][i : i + k + 1]] =                 dkm[data[0][i : i + k + 1]] + 1, dkm[data[1][i : i + k + 1]] + 1
            i = i + 1
    ret = []
    for ky in dkm.keys():
        if thrd < dkm[ky]:
            ret.append(ky)
    return ret

def checkNobranch(adj, v):
    ideg = 0
    odeg = len(adj[v])
    for value in list(adj.values()):
        for u in value:
            if v == u:
                ideg = ideg + 1        
    if odeg == 1 and ideg == 1:
        return True
    return False
    
def cntCont(kms):
    adj = defaultdict(list)
    for km in kms:
        adj[km[:-1]].append(km[1:])
    ret, sts = [], []
    for v in adj:
        if checkNobranch(adj, v) is False:
            sts.append(v)
    i = 0
    while i < len(sts):
        s = sts[i]
        for v in adj[s]:
            nv = v
            p = [s, nv]  
            while True:
                if checkNobranch(adj, nv) is False:
                    break
                nv = adj[nv][0]
                p.append(nv) 
            ret.append(p[0] + ''.join(tp[-1] for tp in p[1:]))
        i = i + 1
    ret = sorted(ret)
    return ret

def genPath(p):
    ret = ""
    for eg in p:
        ret = ret + eg[0]
    ret = ret + p[-1][1:]
    return ret


# In[5]:


if __name__ == "__main__":
    k, fthrd = 36, 2
    datas = parse_reads_file("reads_hw3all_A_3_chr_1.txt")
    if datas is None:
        sys.exit(1)
    cts = sorted(map(genPath, cntCont(preProcess(datas, k, fthrd))))
    with open("r9.txt", 'w') as output_file:
        output_file.write('>reads_hw3all_A_3_chr_1\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(cts))
    #with zipfile.ZipFile(zip_fn, 'w') as myzip:
       # myzip.write("r9.txt")
    


# In[ ]:




