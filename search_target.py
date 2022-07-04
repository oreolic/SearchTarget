#%%
import pandas as pd
from concurrent.futures.process import ProcessPoolExecutor
from datetime import datetime
import os
#%%
class general:        
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
        first = fasta.readline() ## not included first line
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
df =pd.read_csv('Input/input.txt',sep='\t')

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