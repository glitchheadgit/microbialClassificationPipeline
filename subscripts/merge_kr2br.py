#!/usr/bin/env python

# %%
import pandas as pd
import glob
import os
import re
import numpy as np
import bisect
import argparse
# %%
parser = argparse.ArgumentParser(
    prog='bracken_merge',
    description='Merges bracken reports based on kraken_assigned_reads, assembles them with classification from kraken reports',
    epilog=' Made by glitch 	♡( ◡‿◡ )')
parser.add_argument('-i', '--index', default='kraken')
parser.add_argument('-b', '--bracken_dir', required=True)
parser.add_argument('-k', '--kraken_dir', required=True)
parser.add_argument('-n', '--report_name', default='report')
args = parser.parse_args()
idx = args.index
kr_dir=os.path.abspath(args.kraken_dir)
br_dir=os.path.abspath(args.bracken_dir)
fname=args.report_name
# %%
print("Объединяем отчёты bracken...")
b_list = glob.glob(f'{br_dir}/*bracken')
name = re.search(r'([^.]*)\..*', os.path.basename(b_list[0]))[1]
dir = re.search(r'.*\/(.*)\/[^\/]*', b_list[0])[1]
report = pd.read_csv(b_list[0], delimiter='\t')[['name', 'kraken_assigned_reads']]
report = report.rename(columns = {'kraken_assigned_reads':f'{name}_kar_{dir}'})

for file in b_list[1:]:
    name = re.search(r'([^.]*)\..*', os.path.basename(file))[1]
    dir = re.search(r'.*\/(.*)\/[^\/]*', file)[1]
    df_kar = pd.read_csv(file, delimiter='\t')[['name','kraken_assigned_reads']]
    df_kar = df_kar.rename(columns = {'kraken_assigned_reads':f'{name}_kar_{dir}'})
    report = report.merge(df_kar, how="outer", on=["name"])
report = report.rename(columns={'name':'S'})
print("Готово!")

# %%
print("Теперь возьмемся за kraken'a >:3")
columnsn = ['percentage', 'num_fr_cov', 'num_fr_ass', 'rank_code', 'id_num', 'S']
k_list = glob.glob(f'{kr_dir}/*report')
l=len(k_list)
dic = {'G' : [], 'P' : [], 'S' : [], 'K' : [], 'C' : [], 'O' : [], 'F' : [], 'D' : []}
for n, file in enumerate(k_list):
    print(f'[{n+1}/{l}] Обрабатываем {file}...')
    report_k = pd.read_csv(file, delimiter='\t', names=columnsn, header=None, skipinitialspace=True)
    name = re.search(r'([^.]*)\..*', os.path.basename(file))[1]
    for n in report['S']:
        if n in dic['S']:
            continue
        else:
            try:
                ix = min(report_k[report_k['S'].str.endswith(n)].index)
            except:
                continue
            for k in dic.keys():
                dic[k].append(np.nan)
            dic['S'][-1]=n
            bik = bisect.bisect(list(report_k[report_k['rank_code'].str.endswith('D')].index), ix)
            bik = list(report_k[report_k['rank_code'].str.endswith('D')].index)[bik-1]
            dic['D'][-1]=report_k['S'].iloc[bik]
            for key in [k for k in dic.keys() if k not in ['D', 'S']]:
                try:
                    bi = bisect.bisect(list(report_k[report_k['rank_code'].str.endswith(key)].index), ix)
                except:
                    continue
                bi = list(report_k[report_k['rank_code'].str.endswith(key)].index)[bi-1]
                if bi > bik and bi < ix: dic[key][-1]=report_k['S'].iloc[bi]
    print("\033[1A\x1b[2K", end="") # Erase previously printed string

dic = pd.DataFrame.from_dict(dic)

# %%
rm = report.merge(dic, how='left', on='S')

# %%
idx = idx + '_{}'
rm = rm.reset_index(drop=True).rename(idx.format)

# %%
rm_class = rm[['S','G','F','O','C','P','K', 'D']]
rm_kar = rm.drop(['S','G','F','O','C','P','K', 'D'], axis=1)

# %%
rm.to_csv(f'{os.getcwd()}/{fname}.csv')
rm_class.to_csv(f'{os.getcwd()}/{fname}_class.csv')
rm_kar.to_csv(f'{os.getcwd()}/{fname}_kar.csv')
print(f"Отчёты готовы! Они сохранены в текущей директории: {fname}.csv, {fname}_kar.csv, {fname}_class.csv")
