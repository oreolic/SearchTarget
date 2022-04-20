#%%
import pandas as pd
from concurrent.futures.process import ProcessPoolExecutor
from datetime import datetime
import os
#%%
class general:
    def read_fasta(self,num):
        total = ''
        with open('EssentialData/hg38/{}'.format(num) , 'r') as fasta:
            while True:
                a = fasta.readline()
                if not a :
                    break
                else:
                    a  = a.rstrip()
                    total += a

        return total

    def _fasta38(self,num):
        total = ''
        with open('EssentialData/hg38/chr{}.fa'.format(num) , 'r') as fasta:
            a = fasta.readlines()

            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1

        return total

    def _fasta19(self,num):
        total = ''
        with open('EssentialData/hg19/hg19_{}.fasta'.format(num) , 'r') as fasta:
            a = fasta.readlines()
            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1
        return total

    def read_fasta(self,acc):
        hg19 = pd.read_csv('EssentialData/hg19_Accession.txt',sep='\t',index_col=1)
        hg38 = pd.read_csv('EssentialData/hg38_Accession.txt',sep='\t',index_col=1)

        if acc in list(hg19.index):
            num = hg19.loc[acc,'Chromosome'].split()[1]        
            fasta = self._fasta19(num)
            return fasta
        elif acc in list(hg38.index):   
            num = hg38.loc[acc,'Chromosome'].split()[1]           
            fasta = self._fasta38(num)
            return fasta
        else:
            return 'FALSE'
        
    def rt(self,seq):
        dic = {'A':'T','G':'C','T':'A','C':'G','a':'t','g':'c','c':'g','t':'a','-':'-','N':'N'}
        rt = ''
        for i in seq[::-1]:
            rt += dic[i]
        return rt




#%%

def read_fasta(path):
    fa = ''
    with open('EssentialData/hg38/{}'.format(path),'r') as fasta:
        while True:
            a = fasta.readline()
            if not a :
                break
            else:
                a  = a.rstrip()
                fa += a
    return fa

def read_fast19(path):
    fa = ''
    with open('EssentialData/hg19/hg19_{}.fasta'.format(path),'r') as fasta:
        while True:
            a = fasta.readline()
            if not a :
                break
            else:
                a  = a.rstrip()
                fa += a
    return fa


def search_target_fp(fa,tg,index):
    
    if index == -1:
        return 'FALSE'
    else:
        length = len(tg)
        pam = index+length-3

        pre = 14
        post = 28

        reference = fa[pam-length+3-pre:pam+3+post]
        return reference

def search_target_rp(fa,tg,index):
    if index == -1:
        return 'FALSE'
    else:
        pam_rt = index
        pre = 14
        post = 28
        referenc_rt = fa[pam_rt-post:pam_rt+len(tg)+pre]
        reference = general().rt(referenc_rt)
        return reference



# %%
def search_target_in_same_fasta(fa):
    t1 = datetime.now()
    
    df =pd.read_csv('Input/target.txt')
    for idx in df.index:
        tg = df.loc[idx,'Target'].upper()
        df.loc[idx,'Target'] = tg
        if len(tg.split('-')) != 1:
            tg = ''.join(tg.split('-'))
            df.loc[idx,'Target'] = tg

    for idx in df.index:
        tg = df.loc[idx,'Target'].upper()
        ref = search_target_fp(fa,tg)
        if ref == 'FALSE':
            df.loc[idx,'reference'] = ref
        else:
            pass

    result = df[df['Target']!='FALSE']
    t2 = datetime.now()
    print(t2-t1)
    return result



def search_rp_target_in_same_fasta(fa):
    df =pd.read_csv('Input/target.txt')
    for idx in df.index:
        tg = df.loc[idx,'Target'].upper()
        df.loc[idx,'Target'] = tg
        if len(tg.split('-')) != 1:
            tg = ''.join(tg.split('-'))
            df.loc[idx,'Target'] = tg

    for idx in df.index:
        tg = df.loc[idx,'Target'].upper()
        ref = search_target_rp(fa,tg)
        if ref == 'FALSE':
            df.loc[idx,'reference'] = ref
        else:
            pass

    result = df[df['Target']!='FALSE']
    return result

def search_target_in_same_target(tg):
    t1 = datetime.now()
    files = os.listdir('EssentialData/hg38')
    n = 0
    tg_rt = general().rt(tg)
    while n < len(files):
        file = files[n]
        fa = read_fasta(file)

        index = fa.find(tg)
        rt_index = fa.find(tg_rt)

        if index != -1:
            ref = search_target_fp(fa,tg,index)
            return ref
        elif rt_index != -1:
            ref = search_target_rp(fa,tg,rt_index)
            return ref
        else: pass
        n += 1
    

    t2 = datetime.now()
    print(t2-t1)
    return 'FALSE'


            






# %%
files = os.listdir('EssentialData/hg38')
df =pd.read_csv('Input/target.txt',sep='\t')

print('start',datetime.now())
merged = []
with ProcessPoolExecutor(max_workers=40) as executor:
    futs = []
    for idx in df.index:
        tg = df.loc[idx,'Target'].upper()

        tg = tg.rstrip()
        df.loc[idx,'Target'] = tg
        
        if len(tg.split('-')) != 1:
            tg = ''.join(tg.split('-'))
            df.loc[idx,'Target'] = tg
        else:
            pass

        if len(tg.split(',')) != 1:
            tg = tg.split(',')[0]
            df.loc[idx,'Target'] = tg
        general().rt(tg)
        fut  = executor.submit(search_target_in_same_target,tg)
        futs.append((idx,fut))
    
    for tup in futs:
        idx = tup[0]
        ref = tup[1].result()

        df.loc[idx,'reference'] = ref

df.to_csv('Output/find_target.txt',sep='\t',index=None)
print('done',datetime.now())