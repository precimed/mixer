from datetime import datetime
import argparse
import getpass
import glob
import itertools
import ntpath
import numpy as np
import os
import os.path
import pandas as pd
import random
import scipy
import scipy.io as sio
import scipy.linalg
import shutil
import six
import socket
import subprocess
import sys
import sys, traceback
import time
import scipy.stats as st

mysystem = 'MMIL'

# globals (not configurable)

if mysystem == 'MMIL':
    mixer_dir = '/home/oleksandr/github/mixer'
    mixer_lib='/home/oleksandr/github/mixer/src/build/lib/libbgmg.so'
    bfile='<SynGen2>/11015833/bfile_merged/chr'    # for simu_linux
    defvec_hm3 = '<SynGen2>/11015833/bfile_merged/defvec_hapmap3_hardprune_p1.snps'
    ld_file=bfile+'@.run4.ld'
    bim_file=bfile+'@.bim'
    out_prefix_default = '<PROJECTS>/plsa_mixer/simu/run1/simu'
    tmp_prefix_default = os.path.join('/scratch/', 'mixer_simu')
    python_cmd = '~/miniconda3/bin/python3'
    R_cmd = '~/miniconda3/bin/R'
    simu_linux = 'simu_linux'
    load_plink = 'true'

class Params(object):
    def __init__(self, bfile_chr, pi1u, pi2u, pi12, hsq, rg, rep, spow, args):
        self._pi1u = pi1u
        self._pi2u = pi2u
        self._pi12 = pi12
        self._pi1 = pi1u - pi12
        self._pi2 = pi2u - pi12
        self._hsq = str(hsq)
        self._rg = str(rg)
        self._bfile_chr = str(bfile_chr)
        self._rep = str(rep)
        self._spow = spow
        self._args = args

    def num_components(self):
        return str(int(self._pi1 > 0) + int(self._pi2 > 0) + int(self._pi12>0))
    
    def pi_vec(self):
        return [self._pi1, self._pi2, self._pi12]
    
    def trait1_sigsq(self):
        return ' '.join([str(sig) for sig, pi in zip([1, 0, 1], self.pi_vec()) if pi > 0])
    
    def trait2_sigsq(self):
        return ' '.join([str(sig) for sig, pi in zip([0, 1, 1], self.pi_vec()) if pi > 0])
    
    def causal_pi(self):
        return ' '.join(['{:.8e}'.format(pi) for pi in self.pi_vec() if pi > 0])

    def rg(self):
        return ' '.join([str(rgval) for rgval, pi in zip(['0.0', '0.0', self._rg], self.pi_vec()) if pi > 0])

    def traitN_spow(self):
        return ' '.join([str(self._spow) for pi in self.pi_vec() if pi > 0])

    def hsq(self):
        return self._hsq + ' ' + self._hsq

    def out(self, tmp=False):
        return (self._args.tmp_prefix if tmp else self._args.out_prefix ) 

    def plink_out(self, trait, chri):
        return '{out}.pheno.{t}.chr{chri}'.format(chri=chri, t=trait, out=self.out(tmp=True))
  
    def plink_cmd(self, trait, chri):
        extract_snps = '--extract {}'.format(defvec_hm3)
        return load_plink + ' | plink --bfile {bfile}{chri} --assoc {extract} --memory 2048 --allow-no-sex --no-pheno --pheno {out}.simu.pheno --pheno-name {t} --out {o} && gzip -f {o}.qassoc'.format(bfile=bfile, chri=chri, out=self.out(), t=trait, o=self.plink_out(trait, chri), extract=extract_snps)
 
    def cmd(self):
        return  simu_linux + ' --bfile-chr ' + self._bfile_chr + \
               ' --qt --num-traits 2 --hsq ' + self.hsq() + \
               ' --num-components ' + self.num_components() + \
               ' --causal-pi ' + self.causal_pi() + \
               ' --chr-labels ' + ' '.join(self._args.chr2use) + \
               (' --trait1-s-pow ' + self.traitN_spow() if self._spow else '') + \
               (' --trait2-s-pow ' + self.traitN_spow() if self._spow else '') + \
               ' --trait1-sigsq ' + self.trait1_sigsq() + \
               ' --trait2-sigsq ' + self.trait2_sigsq() + \
               ' --rg ' + self.rg() + \
               ' --out ' + self.out() + '.simu'

    def python_mixer_cmd(self, mixer_args):
        return '{} {}/precimed/mixer.py fit \\\n'.format(python_cmd, mixer_dir) + \
        ' \\\n'.join(['\t--{} {}'.format(k, v) for k, v in mixer_args.items()])
    
    def python_bgmg_and_ugmg_common(self):
        mixer_args = {}
        mixer_args['bim-file'] = bim_file
        mixer_args['ld-file'] = ld_file
        mixer_args['chr2use'] = ','.join(self._args.chr2use)
        mixer_args['lib'] = mixer_lib
        mixer_args['threads'] = str(self._args.threads)
        mixer_args['ci-alpha'] = '0.05'
        #mixer_args['fit-sequence'] = 'diffevo-fast'
        #mixer_args['diffevo-fast-repeats'] = '1'
        return mixer_args      

    def python_ugmg_cmd(self, trait):
       mixer_args = self.python_bgmg_and_ugmg_common()
       mixer_args['trait1-file'] = '{}.{}.sumstats.gz'.format(self.out(), trait)
       mixer_args['out'] = '{}.{}'.format(self.out(), trait)
       return self.python_mixer_cmd(mixer_args)

    def python_bgmg_cmd(self):
       mixer_args = self.python_bgmg_and_ugmg_common()
       mixer_args['trait1-file'] = '{}.trait1.sumstats.gz'.format(self.out())
       mixer_args['trait2-file'] = '{}.trait2.sumstats.gz'.format(self.out())
       mixer_args['trait1-params-file'] = '{}.trait1.json'.format(self.out())
       mixer_args['trait2-params-file'] = '{}.trait2.json'.format(self.out())
       mixer_args['out'] = '{}.bgmg'.format(self.out())
       return self.python_mixer_cmd(mixer_args)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = six.moves.reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh, mode):
        self.fh = fh
        self.log_fh = open(fh, mode) if (fh is not None) else None

        # remove error file from previous run if it exists
        try:
            os.remove(fh + '.error')
        except OSError:
            pass

    def system(self, command, dry_run):
        start_time = time.time()
        log.log('>command start time: {T}'.format(T=time.ctime()) )
        self.log(('[DRY-RUN] ' if args.dry_run else '') + command)

        if not dry_run: os.system(command)
 
        time_elapsed = round(time.time()-start_time,2)
        log.log('=command elapsed time: {T}'.format(T=sec_to_str(time_elapsed)))
        log.log('<command end time: {T}'.format(T=time.ctime()) )        

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            self.log_fh.flush()

    def error(self, msg):
        '''
        Print to log file, error file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            with open(self.fh + '.error', 'w') as error_fh:
                error_fh.write(str(msg).rstrip() + '\n')


def parse_args(args):
    parser = argparse.ArgumentParser(description="A helper tool to run simulations with MiXeR.")
    parser.add_argument("--spow", type=str, default=None, help='S parameter in MAF-dependent simulations')
    parser.add_argument("--dry-run", default=False, action="store_true", help='Dry run')
    parser.add_argument("--rep", type=str, default='0', help="simulation repeat")
    parser.add_argument("--h2", type=str, default='0.4', help="heritability")
    parser.add_argument("--rg", type=str, default='0.0', help="correlation of genetic effects")
    parser.add_argument("--pi1u", type=str, default=3e-03, help="polygenicity of the first trait")
    parser.add_argument("--pi2u", type=str, default=3e-03, help="polygenicity of the second trait")
    parser.add_argument("--pifrac", type=str, default=0.4, help="polygenicity of the second trait")
    parser.add_argument("--threads", type=str, default=8, help="number of threads to use in MiXeR computations")
    parser.add_argument("--out-prefix", type=str, default=out_prefix_default)
    parser.add_argument("--tmp-prefix", type=str, default=tmp_prefix_default)
    parser.add_argument("--chr2use", type=str, nargs='+', default=[str(x) for x in range(1, 23)])
    parser.add_argument("--analysis", type=str, default=['simu', 'plink', 'ugmg', 'bgmg'], nargs='*', choices=['simu', 'plink', 'ugmg', 'bgmg'])
    return parser.parse_args(args)


def run_analysis(args, log, param):
    #if os.path.exists(pheno_out):
    #    log.log('skip making {}, file exists'.format(pheno_out))
    #    args.analysis = [x for x in args.analysis if x != 'pheno']

    if 'simu' in args.analysis:
        log.system(param.cmd(), args.dry_run)
        
    if 'plink' in args.analysis:
        task = param.plink_cmd('{1}', '{2}')
        log.system('parallel -j16 "' + task + '" ::: trait1 trait2 ::: {}'.format(' '.join(args.chr2use)), args.dry_run)
        ref = pd.concat([pd.read_csv(bim_file.replace('@', str(chri)), delim_whitespace=True, header=None, names='CHR SNP GP BP A1 A2'.split()) for chri in args.chr2use])
        for trait in ['trait1', 'trait2']:
            plink_out = param.out() + '.{}.sumstats.gz'.format(trait)
            df=pd.concat([pd.read_csv(param.out(tmp=True) + '.pheno.{}.chr{}.qassoc.gz'.format(trait, chri),delim_whitespace=True) for chri in args.chr2use])
            df.reset_index(drop=True, inplace=True)
            df['Z'] = st.norm.ppf(df.P.values*0.5)*np.sign(df.BETA)
            df = pd.merge(df, ref[['SNP', 'A1', 'A2']], on='SNP', how='left')[['SNP', 'A1', 'A2', 'NMISS', 'Z']].rename(columns={'NMISS':'N'})
            df.to_csv(plink_out, index=False, sep='\t')
            log.log('generate {}'.format(plink_out))

    if 'ugmg' in args.analysis:
        log.system(param.python_ugmg_cmd("trait1"), args.dry_run)
        log.system(param.python_ugmg_cmd("trait2"), args.dry_run)

    if 'bgmg' in args.analysis:
        log.system(param.python_bgmg_cmd(), args.dry_run)


def insert_key_to_dictionary_as_list(key, value, df_data):
    if key not in df_data:
        df_data[key] = []
    df_data[key].append(value)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    #if args.seed != None: np.random.seed(args.seed)
   
    param=Params(bfile, float(args.pi1u), float(args.pi2u), np.minimum(float(args.pi2u), float(args.pi1u))*float(args.pifrac), args.h2, args.rg, args.rep, args.spow, args)

    log = Logger(param.out() + '.log', 'a')
    start_time = time.time()

    try:
        defaults = vars(parse_args([]))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = "Call: \n"
        header += './mixer_simu.py \\\n'
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T} by {U}, host {H}'.format(T=time.ctime(), U=getpass.getuser(), H=socket.gethostname()))

        for envvar in 'JOB_ID HOSTNAME QUEUE JOB_NAME JOB_SCRIPT TMP TMPDIR'.split():
            if envvar not in os.environ: continue
            log.log('os.environ["{}"] = {}'.format(envvar, os.environ[envvar]))

        run_analysis(args, log, param)
        if (args.out_prefix != args.tmp_prefix): log.system('rm {}*'.format(param.out(tmp=True)), args.dry_run)

    except Exception:
        log.error( traceback.format_exc() )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
