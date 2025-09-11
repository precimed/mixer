import pandas as pd
import numpy as np

out_prefix = '/home/oleksanf/github/precimed/gsa-mixer/out_example/PGC_SCZ_0518_EUR'
#out_prefix = '/cluster/projects/nn9114k/oleksanf/gsa-mixer/out2/PGC_SCZ_0518_EUR'

# an optional path to gene-sets to be excluded from the results
# a better option is to just filter them from the file passed to --go-file-test in the GSA-MiXeR enrichment model
#prunedoverlap_fname = None
prunedoverlap_fname = '/home/oleksanf/github/precimed/gsa-mixer/reference/gsa-mixer-geneset-annot_10mar2023_prunedoverlap0p8.txt'

def sortcols(df, drop_h2_cols=False, drop_magma_sig_cols=False, drop_mixer_sig_cols=False, sort_col='enrich'):

    if not drop_magma_sig_cols:
        # add MAGMA rank column
        df_sorted = df.sort_values('MAGMA_P', ascending=True).reset_index(drop=True)
        df_sorted['MAGMA_RANK'] = df_sorted.index + 1
        df = df_sorted.copy()

    cols = ['GO', 'GENE', 'NGENES', 'enrich', 'se_enrich', 'MIXER_AIC', 'loglike_diff', 'loglike_df', 'MIXER_P', 'MIXER_Q', 'MIXER_FDR', 'MAGMA_P', 'MAGMA_RANK',  'h2', 'se_h2', 'h2_frac',  'se_h2_frac',  'h2_base_frac', 'se_h2_base_frac', 'excl_GENE', 'excl_enrich', 'se_excl_enrich', 'GENE_LIST']
    df_sorted = df[[c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols]].copy()
    
    if sort_col is not None:
        df_sorted['sort_key'] = -df_sorted[sort_col]
        df_sorted.sort_values(['sort_key'], inplace=True)

    for c in df_sorted.columns:
        if c in 'h2 se_h2 h2_frac se_h2_frac h2_base_frac se_h2_base_frac'.split():
            df_sorted[c] = [('' if pd.isnull(x) else ('{:,.6f}'.format(x))) for x in df_sorted[c].values]
        if c in 'MIXER_AIC loglike_diff'.split():
            df_sorted[c] = [('' if pd.isnull(x) else ('{:,.2f}'.format(x))) for x in df_sorted[c].values]
        if c in 'enrich se_enrich excl_enrich se_excl_enrich'.split():
            df_sorted[c] = [('' if pd.isnull(x) else ('{:,.2f}'.format(x))) for x in df_sorted[c].values]
        if c in 'MAGMA_P MIXER_P MIXER_Q MIXER_FDR'.split():
            df_sorted[c] = [('' if pd.isnull(x) else ('{:,.2e}'.format(x))) for x in df_sorted[c].values]
        if c in ['NGENES']:
            df_sorted[c] = [('' if pd.isnull(x) else int(x)) for x in df_sorted[c].values]
    
    if sort_col is not None: del df_sorted['sort_key']
    
    if drop_magma_sig_cols and ('MAGMA_P' in df_sorted):
        del df_sorted['MAGMA_P']
        
    if drop_mixer_sig_cols:
        del df_sorted['MIXER_AIC']
        del df_sorted['loglike_diff']
        del df_sorted['loglike_df']
    
    if drop_h2_cols:
        del df_sorted['h2']
        del df_sorted['se_h2']
    
    df_sorted=df_sorted.reset_index(drop=True)
    return df_sorted

if __name__ == "__main__":
    print('Reading MAGMA output...')
    df_magma=pd.read_csv(f'{out_prefix}_magma.gsa.out',sep='\s+',comment='#')
    df_magma_gene=pd.read_csv(f'{out_prefix}_magma.step2.genes.out',sep='\s+',comment='#')

    print('Reading GSA-MiXeR output...')
    df = pd.read_csv(f'{out_prefix}_full.go_test_enrich.csv',sep='\t')
    if 'enrich_std' in df.columns: df.rename(columns={'enrich_std':'se_enrich'}, inplace=True)
    if 'h2_std' in df.columns: df.rename(columns={'h2_std':'se_h2'}, inplace=True)

    geneset_pruned = pd.read_csv(prunedoverlap_fname, sep='\t') if (prunedoverlap_fname is not None) else []

    h2_coding_full = df[df['GO']=='coding_genes']['h2'].iloc[0]
    h2_coding_base = df[df['GO']=='coding_genes']['h2_base'].iloc[0]

    df['h2_frac'] = df['h2'] / h2_coding_full
    df['se_h2_frac'] = df['se_h2'] / h2_coding_full
    df['h2_base_frac'] = df['h2_base'] / h2_coding_base
    df.rename(columns={'loglike_aic':'MIXER_AIC'}, inplace=True)
    del df['h2_base']
    del df['snps']

    print('Organize output tables...')
    idx_base_coding = df['GO'].isin(['base', 'coding_genes'])
    idx_geneset_logo = df['GO'].str.contains('_excl_')
    idx_geneset_main = ~idx_geneset_logo & (df['GO'].str.startswith('GOMF_') |
                                            df['GO'].str.startswith('GOBP_') |
                                            df['GO'].str.startswith('GOCC_') |
                                            df['GO'].str.startswith('SYNGO_'))
    idx_geneset_pruned = idx_geneset_main & ~df['GO'].isin(geneset_pruned) & (df['NGENES'] >= 5)
    idx_gene = ~idx_base_coding & ~idx_geneset_logo & ~idx_geneset_main
                
    df_mixer_genes        = df[idx_gene].copy()
    df_mixer_geneset_main = df[idx_geneset_pruned].copy()

    df_mixer_geneset_logo = df[idx_geneset_logo & ~df['enrich'].isnull()].copy()
    df_mixer_geneset_logo.rename(columns={'enrich':'excl_enrich', 'se_enrich':'se_excl_enrich'}, inplace=True)
    df_mixer_geneset_logo['excl_GENE'] = [x.split('_excl_')[1] for x in df_mixer_geneset_logo['GO']]
    df_mixer_geneset_logo['GO'] = [x.split('_excl_')[0] for x in df_mixer_geneset_logo['GO']]
    df_mixer_geneset_logo.reset_index(drop=True, inplace=True)
    df_mixer_geneset_logo = df_mixer_geneset_logo.loc[df_mixer_geneset_logo[['GO', 'excl_enrich', 'se_excl_enrich']].groupby(['GO'])['excl_enrich'].idxmin()]
    df_mixer_geneset_main = pd.merge(df_mixer_geneset_main, df_mixer_geneset_logo[['GO', 'excl_GENE', 'excl_enrich', 'se_excl_enrich']], on=['GO'], how='left')

    df_geneset = pd.merge(df_mixer_geneset_main,
                        df_magma.rename(columns={'FULL_NAME': 'GO', 'P':'MAGMA_P'})[['GO', 'MAGMA_P']],
                        how='left', on=('GO'))

    g2g=pd.DataFrame([x.split('_excl_') for x in df[df['GO'].str.contains('_excl_')]['GO']], columns=['GO', 'GENE'])                        
    df_geneset = pd.merge(df_geneset,
                          g2g.groupby('GO')['GENE'].apply(lambda x:' '.join(sorted(x))).reset_index().rename(columns={'GENE':'GENE_LIST'}),
                          on='GO', how='left')

    df_genes = df_mixer_genes.copy()
    df_genes = pd.merge(df_genes,
                        df_magma_gene.rename(columns={'P':'MAGMA_P', 'GENE':'GO'})[['GO', 'MAGMA_P']],
                        on=['GO'], how='left')
    del df_genes['NGENES']

    if 'GO' in df_genes.columns: df_genes.rename(columns={'GO':'GENE'}, inplace=True)

    with pd.ExcelWriter(f"{out_prefix}_SupplementaryTables.xlsx") as excel_writer:
        print('Output gene-sets, filtered by magma, re-ordered by MiXeR...')
        sortcols(df_geneset[(df_geneset['MAGMA_P'] <  0.05/len(df_magma['FULL_NAME'].unique()))], drop_h2_cols=True, drop_mixer_sig_cols=True, sort_col='enrich'  ).to_excel(excel_writer, sheet_name='ST5', index=False)

        print('Output gene-sets, filtered by MiXeR AIC, ordered by enrich...')
        sortcols(df_geneset[(df_geneset['MIXER_AIC']>0)],
                drop_h2_cols=True, drop_magma_sig_cols=True, sort_col='enrich').to_excel(excel_writer, sheet_name='ST7', index=False)

        print('Output gene-level results...')
        sortcols(df_genes[(df_genes['MIXER_AIC']>0)], 
                drop_magma_sig_cols=True, sort_col='MIXER_AIC').to_excel(excel_writer, sheet_name='ST9', index=False)

        print('Saving output file...')

    print('Done.')
