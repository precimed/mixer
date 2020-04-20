DO_PREPARE_READ_DATA_NOMHC=False
DO_BGMG_DENSITY=False
DO_BGMG_CAUSAL_DENSITY=False
BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY = False
DO_BGMG_CAUSAL_DENSITY_SUPPL=True

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

def load_data(trait1, trait2, censored=False):
    fname1 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait1, trait2)
    fname2 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait2, trait1)
    if os.path.exists(fname2):
        data = json_loads(fname2)
        flip_data = False
    elif os.path.exists(fname1):
        data = json_loads(fname1)
        flip_data = True
    else:
        raise ValueError('missing: {} vs {}'.format(trait1, trait2))
    return data, flip_data

def plot_causal_density(data, flip_data, vmax=1e3, plot_limits=0.025):
    params = data['result']['params']
    sb1, sb2 = params['sig2_beta'][0][2], params['sig2_beta'][1][2]
    pi1, pi2, pi12 = tuple(params['pi_vec'])
    rho =max(min(params['rho_beta'][2], 0.98), -0.98)
    cov = rho * np.sqrt(sb1*sb2)
    factor = 1; sb_null = 1e-7
    #print(sb1, sb2, cov, sb_null, pi1, pi2, pi12)
    rv1 = multivariate_normal([0, 0], [[sb1, 0], [0, sb_null]])
    rv2 = multivariate_normal([0, 0], [[sb_null, 0], [0, sb2]])
    rv12 = multivariate_normal([0, 0], [[sb1, cov], [cov, sb2]])
    grid_step=plot_limits/50
    x, y = np.mgrid[-plot_limits:plot_limits:grid_step, -plot_limits:plot_limits:grid_step]
    pos = np.empty(x.shape + (2,))
    pos[:, :, 0] = x; pos[:, :, 1] = y
    z=factor*1e7*grid_step*grid_step*(pi1*rv1.pdf(pos)+pi2*rv2.pdf(pos)+pi12*rv12.pdf(pos))
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    im=plt.imshow(np.maximum(1,z if flip_data else z.T),interpolation='none', origin='lower', cmap='magma', norm=LogNorm(), vmin=1, vmax=vmax,extent=plot_extent)
    return im

def plot_predicted_zscore(data, flip_data, num_snps):
    density=data['result']['bivariate']['stratified_qq_plot_fit_data']['trait1'][0]
    plot_limits = max(density['pdf_zgrid'])
    plot_step = density['pdf_zgrid'][1]-density['pdf_zgrid'][0]
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    nbins_to_pdfsize_scale = 10 # zDATA histogram is 100x100 grid, while zBGMG histogram is 1000x1000 grid. Therefore counts in zBGMG are 10 times smaller, and we need to adjust for it.
    zBGMG = (nbins_to_pdfsize_scale*nbins_to_pdfsize_scale) * np.array(density['pdf']) * plot_step * plot_step * num_snps
    im=plt.imshow(np.maximum(1, zBGMG.T if flip_data else zBGMG), interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
    plot_limits=15; plt.axis([-plot_limits, plot_limits, -plot_limits, plot_limits])
    return im, zBGMG

if DO_BGMG_DENSITY or DO_BGMG_CAUSAL_DENSITY:
    fig = plt.figure(figsize=(18, 18), dpi=80)

if DO_PREPARE_READ_DATA_NOMHC:
    #read_data_noMHC=None
    if ('read_data_noMHC' not in locals()) or (read_data_noMHC is None):
        read_data_noMHC = {}
        for trait in traits14_ordered:
            fname = format_string_sumstats_noMHC.format(trait)
            print('reading {}...'.format(fname))
            read_data_noMHC[trait] = pd.read_table(fname, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
    
t1vo = 'PGC_SCZ_2014_EUR' # trait1 vs others
if DO_BGMG_DENSITY:
    if ('read_data_noMHC_pairs' not in locals()) or (read_data_noMHC_pairs is None):
        #for trait1, trait2 in list(itertools.combinations(traits4_ordered, 2)):
        read_data_noMHC_pairs = {}
        for trait1, trait2 in [(t1vo, x) for x in traits4_ordered if (x != t1vo)]:
            print('processing {} vs {}'.format(trait1, trait2))
            df1 = read_data_noMHC[trait1]
            df2 = read_data_noMHC[trait2]
            df = merge_z_vs_z(df1, df2)

            if False:  # filter hm3 snps
                hm3_snps = pd.read_table(data_root + r"/MMIL/SUMSTAT/misc/w_hm3.snplist", delim_whitespace=True)
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
        zDATA, _, _ = np.histogram2d(df['Z2'], df['Z1'], bins=100, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
        im=plt.imshow(np.maximum(1,zDATA.T),interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
        if True: plt.xlabel('$z_{'+traits4_short[trait2]+'}$')
        if index==0: plt.ylabel('$z_{'+traits4_short[trait1]+'}$')
        if index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if index!=2: plt.gca().set_visible(False)

if DO_BGMG_CAUSAL_DENSITY:
    for index, (trait1, trait2) in enumerate([(t1vo, x) for x in traits4_ordered if (x != t1vo)]):
        ti1 = trait4_index_map[trait1]
        ti2 = trait4_index_map[trait2]
        num_snps_data = len(read_data_noMHC_pairs[(trait1, trait2)])
        plt.subplot(3,3,7 + index)

        data, flip_data = load_data(trait1, trait2, censored=False)
        im = plot_causal_density(data, flip_data)
        if True: plt.xlabel('$\\beta_{'+traits4_short[trait2]+'}$')
        if index==0: plt.ylabel('$\\beta_{'+traits4_short[trait1]+'}$')
        if index!=0: plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
        if index!=2: plt.gca().set_visible(False)

        # BGMG-estimated density of Z scores
        plt.subplot(3,3,4 + index)
        im, zBGMG = plot_predicted_zscore(data, flip_data, num_snps_data)
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
    savefig(figures_folder, '{}'.format(figure_name))

if BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY:
    #for file in glob.glob(r'H:\work\run_simu_bgmg_paper_examples\*.pheno'):
    fig = plt.figure(figsize=(18, 12), dpi=80)
    phenos = [folder_simu_bgmg_paper_examples + r'/simu_h2=0.4_rg=0.0_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=0.0000e+00_rep=1_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.pheno',  
              folder_simu_bgmg_paper_examples + r'/simu_h2=0.4_rg=0.99_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=5.0000e-05_rep=1_tag1=customPolygenicOverlapAt0p5_tag2=evenPolygenicity.pheno', 
              folder_simu_bgmg_paper_examples + r'/simu_h2=0.4_rg=0.0_pi1u=1.0000e-04_pi2u=1.0000e-04_pi12=5.0000e-05_rep=1_tag1=customPolygenicOverlapAt0p5_tag2=evenPolygenicity.pheno']
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

    savefig(figures_folder, 'BGMG_SIMU_DENSITY_WITH_CAUSAL_DENSITY')
    
if DO_BGMG_CAUSAL_DENSITY_SUPPL:
    
    causal_density_suppl_tasks = []
    for trait1 in traits_ordered:
        causal_density_suppl_tasks.append((traits14_short[trait1], 1e3, 0.025, ([(trait1, trait2) for trait2 in traits_ordered if trait2 != trait1])))
    causal_density_suppl_tasks.append(('immuno', 1e2, 0.1, list(itertools.combinations(traitsIM_ordered, 2))))
    causal_density_suppl_tasks.append(('antro', 1e2, 0.025, list(itertools.combinations(traitsAM_ordered, 2))))
    censored=False
    
    for task_name, vmax, plot_limits_causal, causal_density_suppl_sub_tasks in causal_density_suppl_tasks:
        fig=plt.figure(figsize=(24, 42), dpi=80)
        for index, (trait1_orig, trait2_orig) in enumerate(causal_density_suppl_sub_tasks):
            trait1 = trait1_orig; trait2=trait2_orig
            print('{} vs {}...'.format(trait1_orig, trait2_orig))
            plt.subplot(len(causal_density_suppl_sub_tasks), 3, 3*index + 3)

            data, flip_data = load_data(trait1, trait2, censored)
            im = plot_causal_density(data, flip_data, vmax, plot_limits_causal)
            if True: plt.xlabel('$\\beta_{'+traits14_short[trait2]+'}$')
            if True: plt.ylabel('$\\beta_{'+traits14_short[trait1]+'}$')
            #if index!=0: plt.gca().get_yaxis().set_visible(False)
            plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
            #if index!=2: plt.gca().set_visible(False)

            df = merge_z_vs_z(read_data_noMHC[trait1], read_data_noMHC[trait2])
            plt.subplot(len(causal_density_suppl_sub_tasks), 3, 3*index + 1)
            plot_limits=15; plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
            z, _, _ = np.histogram2d(df['Z2'], df['Z1'], bins=100, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
            im=plt.imshow(np.maximum(1,z.T),interpolation='none', origin='lower', cmap='hot', norm=LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
            plt.xlabel('$z_{'+traits14_short[trait2]+'}$')
            plt.ylabel('$z_{'+traits14_short[trait1]+'}$')
            plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

            # BGMG-estimated density of Z scores
            plt.subplot(len(causal_density_suppl_sub_tasks), 3, 3*index + 2)
            im, zBGMG = plot_predicted_zscore(data, flip_data, len(df))
            if True: plt.xlabel('$\\hat z_{'+traits14_short[trait2]+'}$')
            if True: plt.ylabel('$\\hat z_{'+traits14_short[trait1]+'}$')
            #if index!=0: plt.gca().get_yaxis().set_visible(False)
            plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))
            #if index!=2: plt.gca().set_visible(False)

        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        savefig(figures_folder, 'BGMG_CAUSAL_DENSITY_SUPPL_{}'.format(task_name) + (' (censored)' if censored else '') )

EVALUATE_HAPGEN = False
if EVALUATE_HAPGEN:
    if 0:
        df_joint=pd.read_table('/home/oleksanf/vmshare/data/bfile_merged_relcheck/1kG(x)_and_hapgen(y)_r2.csv',sep='\t')
        df=pd.read_table('/home/oleksanf/vmshare/data/bfile_merged_relcheck/1kG(x)_and_hapgen(y)_ld.csv', delim_whitespace=True)
        df_pca=pd.read_table('/home/oleksanf/vmshare/data/bfile_merged_relcheck/10k_indep_pca.eigenvec',header=None, delim_whitespace=True,names=['id1', 'id2'] + ['pca{}'.format(i) for i in range(1, 21)])

    fig = plt.figure(figsize=(24, 24), dpi=80)
    plt.subplot(2,2,1)
    n=plt.hist([df['TLD_x'].values, df['TLD_y'].values*nsnps_LDSR/snps_HAPGEN], np.linspace(0, 25, 25), label=['1kG EUR', 'HapGen'])
    plt.legend(loc='upper right')

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('LD score')
    plt.ylabel('#SNPs')
    plt.title('A', loc='right')

    plt.subplot(2,2,2)
    plt.hist([df['TLD_x'].values, df['TLD_y'].values*nsnps_LDSR/snps_HAPGEN], np.linspace(0, 500, 25), label=['1kG EUR', 'HapGen'])
    plt.legend(loc='upper right')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('LD score')
    plt.ylabel('#SNPs')
    plt.title('B', loc='right')

    d1=df_joint[df_joint['R2_x'].isnull()]
    d2=df_joint[df_joint['R2_y'].isnull()]
    d3=df_joint[(~df_joint['R2_y'].isnull()) &(~df_joint['R2_x'].isnull())]
    plt.subplot(2,2,3)
    plt.scatter(d1.SNP.values, d1.TAG.values, s=2, c='r')
    plt.scatter(d2.SNP.values, d2.TAG.values, s=1, c='silver')
    plt.scatter(d3.SNP.values, d3.TAG.values, s=1, c='k')
    plt.ylabel('SNP index')
    plt.xlabel('SNP index')
    plt.title('C', loc='right')

    plt.subplot(2,2,4)
    plt.scatter(df_pca['pca1'].values, df_pca['pca2'].values,s=1 )
    plt.ylabel('PCA component #1')
    plt.xlabel('PCA component #2')
    plt.title('D', loc='right')

    savefig(figures_folder, 'evaluate_hapgen', ['png'])

if 0: # data prep
    import precimed
    import precimed.mixer
    libbgmg = precimed.mixer.LibBgmg('/home/oleksanf/github/mixer/src/build/lib/libbgmg.so', context_id=1)
    libbgmg.init_log("/home/oleksanf/github/mixer/testlog5.log")
    libbgmg.log_message('Test log message succeeded?')
    libbgmg.dispose()

    bim_file = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim'
    frq_file = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink_freq/1000G.EUR.QC.@.frq'
    plink_ld_bin = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.p05_SNPwind50k.ld.bin'
    chr_labels = [1] # list(range(1, 23))
    trait1_file = ''; #'/home/oleksanf/vmshare/data/MMIL/SUMSTAT/TMP/ldsr/PGC_SCZ_2014_EUR.sumstats.gz'
    exclude = ''; extract = ''
    libbgmg.init(bim_file, frq_file, chr_labels, trait1_file, '', exclude, extract);
    print(libbgmg)

    options=[('r2min', 0.05), ('kmax', 100), ('max_causals', 0.03*libbgmg.num_snp), ('num_components', 1), 
             ('cache_tag_r2sum', False), ('threads', 6), ('seed', None), ('z1max', None)]
    for opt, val in options: libbgmg.set_option(opt, val)

    for chr_label in chr_labels: 
        libbgmg.set_ld_r2_coo_from_file(plink_ld_bin.replace('@', str(chr_label)))
        libbgmg.set_ld_r2_csr(chr_label);

    libbgmg.set_option('diag', 0)
    libbgmg_hapgen = precimed.mixer.LibBgmg('/home/oleksanf/github/mixer/src/build/lib/libbgmg.so', context_id=2)
    libbgmg_hapgen.dispose()

    bim_file = '/home/oleksanf/vmshare/data/bfile_merged/chr@.bim'
    frq_file = '/home/oleksanf/vmshare/data/bfile_merged/chr@.frq'
    plink_ld_bin = '/home/oleksanf/vmshare/data/hapgen_ldmat2_plink/bfile_merged_ldmat_p01_SNPwind50k_chr@.ld.bin'
    chr_labels = [1] # list(range(1, 23))
    trait1_file = ''
    exclude = ''; extract = ''
    libbgmg_hapgen.init(bim_file, frq_file, chr_labels, trait1_file, '', exclude, extract);
    print(libbgmg_hapgen)

    options=[('r2min', 0.05), ('kmax', 100), ('max_causals', 0.03*libbgmg_hapgen.num_snp), ('num_components', 1), 
             ('cache_tag_r2sum', False), ('threads', 6), ('seed', None), ('z1max', None)]
    for opt, val in options: libbgmg_hapgen.set_option(opt, val)

    for chr_label in chr_labels: 
        libbgmg_hapgen.set_ld_r2_coo_from_file(plink_ld_bin.replace('@', str(chr_label)))
        libbgmg_hapgen.set_ld_r2_csr(chr_label);

    libbgmg_hapgen.set_option('diag', 0)
    
    import pandas as pd

    chrlist = [1]

    df_hapgen = pd.concat([pd.read_table('/home/oleksanf/vmshare/data/bfile_merged/chr{}.bim'.format(chri),delim_whitespace=True,header=None, names='CHR SNP GP BP A1 A2'.split()) for chri in chrlist])
    df_hapgen.reset_index(drop=True, inplace=True)
    df_hapgen['index']=df_hapgen.index
    df_hapgen['TLD'] = libbgmg_hapgen.ld_tag_r2_sum

    df = pd.concat([pd.read_table('/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim'.format(chri),delim_whitespace=True,header=None, names='CHR SNP GP BP A1 A2'.split()) for chri in chrlist])
    df.reset_index(drop=True, inplace=True)
    df['index']=df.index
    df['TLD'] = libbgmg.ld_tag_r2_sum

    df_merge = pd.merge(df, df_hapgen, on='SNP', how='inner')

    df.to_csv('/home/oleksanf/vmshare/data/bfile_merged_relcheck/1kG_EUR_ld.csv',sep='\t',index=False)
    df_hapgen.to_csv('/home/oleksanf/vmshare/data/bfile_merged_relcheck/hapgen_ld.csv',sep='\t',index=False)
    df_merge.to_csv('/home/oleksanf/vmshare/data/bfile_merged_relcheck/1kG(x)_and_hapgen(y)_ld.csv',sep='\t',index=False)
    
    index_y = df_merge[(df_merge.CHR_x==1) & (df_merge.BP_x<2e6)]['index_y'].values
    index_x = df_merge[(df_merge.CHR_x==1) & (df_merge.BP_x<2e6)]['index_x'].values

    [hSNP, hTAG, hR2] = libbgmg_hapgen.get_ld_r2_chr(chr_label=1)
    [kgSNP, kgTAG, kgR2] = libbgmg.get_ld_r2_chr(chr_label=1)
    df_hR2 = pd.DataFrame({'SNP':hSNP, 'TAG':hTAG, 'R2':hR2})
    df_kgR2 = pd.DataFrame({'SNP':kgSNP, 'TAG':kgTAG, 'R2':kgR2})

    r2thresh = 0.2

    df_kgR2=df_kgR2[(df_kgR2.SNP <= max(index_x)) & (df_kgR2.TAG <= max(index_x)) & (df_kgR2.R2 >= r2thresh)]
    df_hR2=df_hR2[(df_hR2.SNP <= max(index_y)) & (df_hR2.TAG <= max(index_y)) & (df_hR2.R2 >= r2thresh)]

    df_hR2 = df_hR2[df_hR2.SNP.isin(index_y) & df_hR2.TAG.isin(index_y)]
    df_kgR2 = df_kgR2[df_kgR2.SNP.isin(index_x) & df_kgR2.TAG.isin(index_x)]

    for col in ['SNP', 'TAG']:
        df_kgR2[col] = df_kgR2[col].map({v:i for (i,v) in enumerate(index_x)})
        df_hR2[col] = df_hR2[col].map({v:i for (i,v) in enumerate(index_y)})

    df_joint = pd.merge(df_kgR2, df_hR2,how='outer',on=['SNP', 'TAG'])
    df_joint.to_csv('/home/oleksanf/vmshare/data/bfile_merged_relcheck/1kG(x)_and_hapgen(y)_r2.csv',sep='\t',index=False)