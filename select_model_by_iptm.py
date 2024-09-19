import json
import glob
import pandas as pd
from functools import reduce
import sys
from ast import literal_eval
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches

"""
v4: First parameter changed from dir to list
"""

#in_dir = sys.argv[1]
in_file_list = sys.argv[1]
out_prefix = sys.argv[2]

to_check = ['alpha', 'beta', 'mhc', 'ep']

def json2table(summary_file):
    with open(summary_file) as F:
        prefix1 = summary_file.split('/')[-2]
        prefix2 = summary_file.split('_')[-1].replace('.json', '')
        prefix = prefix1 + '_' + prefix2
        ##
        #pp = prefix1.split("_")
        #if pp[2] == "":
            #sample = "_".join(pp[0:2]) + "pos"
        #else:
        #sample = "_".join(prefix1.split("_")[1:3])
        sample = prefix1
        print(sample)
        summary_js = json.load(F)
        summary_df = pd.DataFrame.from_dict(summary_js)
        summary_df['model'] = prefix
        summary_df['model_chain'] = to_check    #this is the reason why order check is required.
        summary_df['sample'] = sample
    return summary_df

def file2list(tfile):
    res = []
    with open(tfile) as F:
        for line in F:
            res.append(line.strip())
    return res

#in_pat = f'{in_dir}/*/*summary_confidences_*.json'
#summary_files = glob.glob(in_pat)
summary_files = file2list(in_file_list)
print(f'Found {len(summary_files)} summary files')

df_list = []
for f in summary_files:
    print(f)
    dt = json2table(f)
    df_list.append(dt)

print('\nmerge ...')
df = reduce(lambda x,y: pd.concat([x,y]), df_list)
df.to_csv(out_prefix+'.xls', sep='\t', index=False)


print('\nmetrics: ep-alpha + ep-beta ...')
df_ep = df[df['model_chain'] == 'ep']    #get ep chain
#df_ep['chain_pair_iptm'] = df_ep['chain_pair_iptm'].apply(lambda x: literal_eval(x))
df_ep['iptm__ep_a_b'] = df_ep['chain_pair_iptm'].apply(lambda x: x[0]+x[1])
df_ep.to_csv(out_prefix+'.ep.xls', sep='\t', index=False)

def add_grp(x):
    xs = x.split('_')
    if len(xs) == 2:
        return 'pos'
    elif xs[2] == "":
        return 'pos'
    else:
        return xs[2]

print('\nboxplot ...')
df_ep['p_ab_iptm'] = df_ep['chain_pair_iptm'].apply(lambda x: x[0:2])
df_ep['p_ab'] = [['alpha', 'beta']] * len(df_ep)

df_ep['grp1'] = df_ep['sample'].apply(lambda x: add_grp(x))
df_ep['grp2'] = df_ep['sample'].apply(lambda x: "_".join(x.split("_")[1:]))

#breakpoint()

dt = df_ep[['p_ab_iptm', 'grp2']].set_index('grp2')
dt = pd.DataFrame(dt['p_ab_iptm'].to_list(), columns=['a', 'b'], index=dt.index)
dt_m = dt.reset_index().groupby('grp2').median().reset_index()

def add_label(x):
    xx = x.split("_")
    if len(xx) == 1:
        label = '1'
    else:
        label = '0'
    return label

dt_m['label'] = dt_m['grp2'].apply(lambda x: add_label(x))
dt_m.to_csv(out_prefix+'.median_label.xls', sep='\t', index=False)

'''hightlight every 6 ticks
##########################################################################
# 20240823

## split too large dataframe
human = file2list('/work/renshuaibing/db/TCR3d/mhc_i_241/202110_202405/human20/human20_rm8SHI.list')
human.append('8SHI')
hu = [ i.lower() for i in human ]

dt = df_ep.explode(['p_ab_iptm', 'p_ab']).sort_values("grp2")
dt1=dt[dt['grp2'].str.contains('|'.join(hu))]
dt2=dt[~dt['grp2'].str.contains('|'.join(hu))]

## boxplot hightlight, mannually replace dt1 or dt2
plt.figure(figsize=(20,10))
sns.boxplot(x="grp2", y="p_ab_iptm", hue="p_ab", data=dt2, fliersize=1)

ax = plt.gca()
xticks = ax.get_xticks()
highlight_indices = [i for i in range(len(xticks)) if i % 6 == 0]

for index in highlight_indices:
    rect = patches.Rectangle((xticks[index] - 0.5, ax.get_ylim()[0]), 1, ax.get_ylim()[1] - ax.get_ylim()[0],
                             linewidth=0, edgecolor=None, facecolor='yellow', alpha=0.3, zorder=0)
    ax.add_patch(rect)

plt.axhline(y=0.8, color='blue', linestyle='--')
plt.xticks(rotation=90, fontsize=8)

#breakpoint()
plt.savefig(out_prefix+'.ep_2.png', bbox_inches='tight')
##########################################################################
'''

'''
plt.figure(figsize=(16,8))
#sns.boxplot(x="sample", y="p_ab_iptm", hue="p_ab", data=df_ep.explode(['p_ab_iptm', 'p_ab']).sort_values("sample"))
#sns.boxplot(x="grp", y="p_ab_iptm", hue="p_ab", data=df_ep.explode(['p_ab_iptm', 'p_ab']).sort_values("sample"))
sns.boxplot(x="grp1", y="p_ab_iptm", hue="p_ab", data=df_ep.explode(['p_ab_iptm', 'p_ab']).sort_values("grp1"))
plt.axhline(y=0.8)
plt.xticks(rotation=90)

#ax = plt.gca()
#rect = patches.Rectangle((22-0.5, ax.get_ylim()[0]), 1, ax.get_ylim()[1]-ax.get_ylim()[0], linewidth=0, edgecolor=None, facecolor='yellow', alpha=0.3, zorder=0)
#ax.add_patch(rect)

plt.savefig(out_prefix+'.ep.png', bbox_inches='tight')
'''
