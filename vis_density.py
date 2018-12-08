DO_BGMG_DENSITY=False
DO_BGMG_CAUSAL_DENSITY=True
BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY = False

def merge_z_vs_z(df1, df2):
    _N_CHR = 22
    # complementary bases
    COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # bases
    BASES = COMPLEMENT.keys()
    # true iff strand ambiguous
    STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                        for x in itertools.product(BASES, BASES)
                        if x[0] != x[1]}
    # SNPS we want to keep (pairs of alleles)
    VALID_SNPS = {x for x in map(lambda y: ''.join(y), itertools.product(BASES, BASES))
                  if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
    # T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
    MATCH_ALLELES = {x for x in map(lambda y: ''.join(y), itertools.product(VALID_SNPS, VALID_SNPS))
                     # strand and ref match
                     if ((x[0] == x[2]) and (x[1] == x[3])) or
                     # ref match, strand flip
                     ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                     # ref flip, strand match
                     ((x[0] == x[3]) and (x[1] == x[2])) or
                     ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
    # T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
    FLIP_ALLELES = {''.join(x):
                    ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                    # strand flip
                    ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                    for x in MATCH_ALLELES} 

    df1 = df1[['SNP', 'A1', 'A2', 'Z']].rename(columns={'Z': 'Z1'}).copy()
    df2 = df2[['SNP', 'A1', 'A2', 'Z']].rename(columns={'Z': 'Z2', 'A1': 'A1x', 'A2': 'A2x'}).copy()
    df = pd.merge(df1, df2, how='inner', on='SNP')
    df = df.dropna(how='any')
    alleles = df.A1 + df.A2 + df.A1x + df.A2x
    df = df[alleles.apply(lambda y: y in MATCH_ALLELES)]
    alleles = df.A1 + df.A2 + df.A1x + df.A2x
    flip_status = alleles.apply(lambda y: FLIP_ALLELES[y])
    df.Z2 *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    df = df.drop(['A1', 'A1x', 'A2', 'A2x'], axis=1)
    return df

def load_data(trait1, trait2):
    fname1 = format_string_BGMG_fit.format(trait1, trait2)
    fname2 = format_string_BGMG_fit.format(trait2, trait1)
    if os.path.exists(fname2):
        data = json.loads(open(fname2).read())
        flip_data = False
    elif os.path.exists(fname1):
        data = json.loads(open(fname1).read())        
        flip_data = True
    else:
        raise ValueError('missing: {} vs {}'.format(trait1, trait2))
    return data, flip_data

def plot_causal_density(data, flip_data):
    params = data['result']['bivariate']['params']
    sb1, sb2 = params['sig2_beta'][0][2], params['sig2_beta'][1][2]
    pi1, pi2, pi12 = tuple(params['pi_vec'])
    rho =max(min(params['rho_beta'][2], 0.98), -0.98)
    cov = rho * np.sqrt(sb1*sb2)
    factor = 1; sb_null = 1e-7
    rv1 = multivariate_normal([0, 0], [[sb1, 0], [0, sb_null]])
    rv2 = multivariate_normal([0, 0], [[sb_null, 0], [0, sb2]])
    rv12 = multivariate_normal([0, 0], [[sb1, cov], [cov, sb2]])
    plot_limits=0.025
    grid_step=plot_limits/50
    x, y = np.mgrid[-plot_limits:plot_limits:grid_step, -plot_limits:plot_limits:grid_step]
    pos = np.empty(x.shape + (2,))
    pos[:, :, 0] = x; pos[:, :, 1] = y
    z=factor*1e7*grid_step*grid_step*(pi1*rv1.pdf(pos)+pi2*rv2.pdf(pos)+pi12*rv12.pdf(pos))
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    im=plt.imshow(np.maximum(1,z if flip_data else z.T),interpolation='none', origin='lower', cmap='magma', norm=LogNorm(), vmin=1, vmax=1e3,extent=plot_extent)
    return im

def plot_predicted_zscore(data, flip_data):
    density=data['result']['bivariate']['stratified_qq_plot_fit_data']['trait1'][0]
    plot_limits = max(density['pdf_zgrid'])
    plot_step = density['pdf_zgrid'][1]-density['pdf_zgrid'][0]
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    z = np.array(density['pdf']) * plot_step * plot_step * 1e7
    im=plt.imshow(np.maximum(1, z.T if flip_data else z), interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
    plot_limits=15; plt.axis([-plot_limits, plot_limits, -plot_limits, plot_limits])
    return im

if DO_BGMG_DENSITY or DO_BGMG_CAUSAL_DENSITY:
    fig = plt.figure(figsize=(18, 18), dpi=80)

t1vo = 'PGC_SCZ_2014_EUR_qc' # trait1 vs others
if DO_BGMG_DENSITY:
    if ('read_data_noMHC' not in locals()) or (read_data_noMHC_pairs is None):
        read_data_noMHC = {}
        for trait in traits_ordered:
            fname = folder_SUMSTAT + r'\STD\{}_noMHC.csv.gz'.format(trait)
            print('reading {}...'.format(fname))
            read_data_noMHC[trait] = pd.read_table(fname, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])

    if ('read_data_noMHC_pairs' not in locals()) or (read_data_noMHC_pairs is None):
        #for trait1, trait2 in list(itertools.combinations(traits4_ordered, 2)):
        for trait1, trait2 in [(t1vo, x) for x in traits4_ordered if (x != t1vo)]:
            print('processing {} vs {}'.format(trait1, trait2))
            df1 = read_data_noMHC[trait1]
            df2 = read_data_noMHC[trait2]
            df = merge_z_vs_z(df1, df2)

            if False:  # filter hm3 snps
                hm3_snps = pd.read_table(folder_SUMSTAT + r"\LDSR\w_hm3.snplist", delim_whitespace=True)
                df = pd.merge(df, hm3_snps, how='inner', on='SNP')
                use_hm3_snps = True
            else:
                use_hm3_snps = False

            read_data_noMHC_pairs[(trait1, trait2)] = df

        for trait1, trait2 in list(read_data_noMHC_pairs.keys()):
            read_data_noMHC_pairs[(trait2, trait1)] = read_data_noMHC_pairs[(trait1, trait2)].copy().rename(columns={'Z1': 'Z2', 'Z2':'Z1'})
    
    for index, (trait1, trait2) in enumerate([(t1vo, x) for x in traits4_ordered if (x != t1vo)]):
        # density of Z scores
        ti1 = trait4_index_map[trait1]
        ti2 = trait4_index_map[trait2]
        plot_limits=15
        plt.subplot(3,3,1 + index)
        df = read_data_noMHC_pairs[(trait1, trait2)]
        plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
        z, _, _ = np.histogram2d(df['Z2'], df['Z1'], bins=100, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
        im=plt.imshow(np.maximum(1,z.T),interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
        if True: plt.xlabel('$z_{'+traits4_short[trait2]+'}$')
        if index==0: plt.ylabel('$z_{'+traits4_short[trait1]+'}$')
        if index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if index!=2: plt.gca().set_visible(False)

if DO_BGMG_CAUSAL_DENSITY:
    for index, (trait1, trait2) in enumerate([(t1vo, x) for x in traits4_ordered if (x != t1vo)]):
        ti1 = trait4_index_map[trait1]
        ti2 = trait4_index_map[trait2]
        plt.subplot(3,3,7 + index)

        data, flip_data = load_data(trait1, trait2)
        im = plot_causal_density(data, flip_data)
        if True: plt.xlabel('$\\beta_{'+traits4_short[trait2]+'}$')
        if index==0: plt.ylabel('$\\beta_{'+traits4_short[trait1]+'}$')
        if index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if index!=2: plt.gca().set_visible(False)

        # BGMG-estimated density of Z scores
        plt.subplot(3,3,4 + index)
        im = plot_predicted_zscore(data, flip_data)
        plot_limits=15; plt.axis([-plot_limits, plot_limits, -plot_limits, plot_limits])
        if True: plt.xlabel('$\\hat z_{'+traits4_short[trait2]+'}$')
        if index==0: plt.ylabel('$\\hat z_{'+traits4_short[trait1]+'}$')
        if index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if index!=2: plt.gca().set_visible(False)

if DO_BGMG_DENSITY and DO_BGMG_CAUSAL_DENSITY:
    figure_name = 'BGMG_DENSITY{}_WITH_CAUSAL_DENSITY'.format('_hm3' if use_hm3_snps else '')
elif DO_BGMG_DENSITY:
    figure_name = 'BGMG_DENSITY{}'.format('_hm3' if use_hm3_snps else '')
else:
    figure_name = 'BGMG_CAUSAL_DENSITY'

if DO_BGMG_DENSITY or DO_BGMG_CAUSAL_DENSITY:
    fig.subplots_adjust(hspace=0.20, wspace=0.20)
    #plt.tight_layout()
    savefig(os.path.join(figures_folder, '{}.svg'.format(figure_name)))
    savefig(os.path.join(figures_folder, '{}.png'.format(figure_name)))

if BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY:
    #for file in glob.glob(r'H:\work\run_simu_bgmg_paper_examples\*.pheno'):
    fig = plt.figure(figsize=(18, 12), dpi=80)
    phenos = [folder_simu_bgmg_paper_examples + r'\simu_h2=0.4_rg=0.0_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=0.0000e+00_rep=1_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.pheno',  
              folder_simu_bgmg_paper_examples + r'\simu_h2=0.4_rg=0.99_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=5.0000e-05_rep=1_tag1=customPolygenicOverlapAt0p5_tag2=evenPolygenicity.pheno', 
              folder_simu_bgmg_paper_examples + r'\simu_h2=0.4_rg=0.0_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=5.0000e-05_rep=1_tag1=customPolygenicOverlapAt0p5_tag2=evenPolygenicity.pheno']
    for file_index, file in enumerate(phenos):
        print('loading {}'.format(file))
        file = os.path.splitext(file)[0]
        fname = os.path.basename(file)
        causals = file + '.{}.causals'
        qassoc = file + '.pheno.trait{}.chr{}.qassoc.gz'

        df_causals = [pd.read_table(causals.format(trait_index), sep='\t') for trait_index in range(1, 3)]
        df_qassoc = [pd.concat([pd.read_table(qassoc.format(trait_index, chr_index), delim_whitespace=True) for chr_index in range(1, 23)]) for trait_index in range(1, 3)]

        if (df_causals[1].SNP!=df_causals[0].SNP).any(): raise("SNPs are different - need pd.merge")

        for trait_index in range(0,2):
            BETA_cols = [c for c in df_causals[trait_index] if c.startswith('BETA_')]
            if np.max(np.array([[df_causals[trait_index][c] != 0] for c in BETA_cols]).sum(axis=0)) > 1: raise('error')
            df_causals[trait_index]['BETA'] = np.squeeze(np.array([[df_causals[trait_index][c]] for c in BETA_cols]).sum(axis=0))

        for trait_index in range(0, 2):
            df_qassoc[trait_index]['Z'] = np.divide(df_qassoc[trait_index]['BETA'].values, df_qassoc[trait_index]['SE'].values)

        # plot z-score vs z-score
        plt.subplot(2,3,4+file_index)
        plot_limits=15
        plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
        z, _, _ = np.histogram2d(df_qassoc[0].Z, df_qassoc[1].Z, bins=100, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
        im=plt.imshow(1+z.T,interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e5,extent=plot_extent)
        plt.axis(plot_extent)
        plt.xlabel('$z_{trait1}$')
        if file_index==0: plt.ylabel('$z_{trait2}$')
        if file_index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if file_index!=2: plt.gca().set_visible(False)

        # plot beta vs beta
        plt.subplot(2,3,1+file_index)
        plot_limits=0.15
        plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits];
        sb_null = 1e-7; x1, y1 = np.random.multivariate_normal([0, 0], [[sb_null, 0], [0, sb_null]], len(df_causals[0])).T
        z, _, _ = np.histogram2d(df_causals[0]['BETA'].values+x1, df_causals[1]['BETA'].values+y1, bins=100, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
        im=plt.imshow(1+z.T,interpolation='none', origin='lower', cmap='magma', norm=LogNorm(), vmin=1, vmax=30,extent=plot_extent)
        plt.axis(plot_extent)
        plt.title(['(A) Independent traits', '(B) Polygenic overlap\nwith genetic correlation', '(C) Polygenic overlap\nw/o genetic correlation'][file_index])
        plt.xlabel('$\\beta_{trait1}$')
        if file_index==0: plt.ylabel('$\\beta_{trait2}$')
        if file_index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if file_index!=2: plt.gca().set_visible(False)

    fig.subplots_adjust(hspace=0.15, wspace=0.15)
    #plt.tight_layout()

    #savefig(os.path.join(figures_folder, '{}.svg'.format(fname)))
    #savefig(os.path.join(figures_folder, '{}.png'.format(fname)))
    savefig(os.path.join(figures_folder, '{}.svg'.format('BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY')))
    savefig(os.path.join(figures_folder, '{}.png'.format('BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY')))
