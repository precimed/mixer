import pandas as pd
import numpy as np
import seaborn as sns
fname='/home/oleksanf/vmshare/data/mixer_test/perf_cluster7/PGC_SCZ_2014_EUR.csv'
#fname='/home/oleksanf/vmshare/data/mixer_test/perf_cluster7/PGC_SCZ_2014_EUR_vs_SSGAC_EDU_2018_no23andMe.csv'
#fname='/home/oleksanf/vmshare/data/mixer_test/perf/PGC_SCZ_2014_EUR.csv'
#fname='/home/oleksanf/vmshare/data/mixer_test/perf/PGC_SCZ_2014_EUR_vs_SSGAC_EDU_2018_no23andMe.csv'
df=pd.read_csv(fname, sep='\t')
df['scale'] = np.nan
for index, row in df.iterrows():
    scale_index = ((df['kmax']==row.kmax) & (df['costcalc']==row.costcalc) & (df['threads']==1))
    df.loc[index, 'scale'] = df[scale_index]['time_sec'].values[0]
df['scale'] = np.divide(df['scale'].values, df['time_sec'].values)

dfc=df.copy()
dfc['costcalc'][:] ='true'
dfc['scale']=dfc['threads']
dfc['time_sec'][:]=np.nan
dfc['cost'][:]=np.nan
dfc.drop_duplicates(['threads', 'kmax', 'costcalc'], inplace=True)
df=pd.concat([dfc, df])

sns.catplot(x='threads', y='scale', hue='costcalc', col='kmax', data=df, kind='bar')
