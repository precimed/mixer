import os.path
import sys
import glob

design = 'run1'

if design == 'run1':
    out_prefix = '<PROJECTS>/plsa_mixer/simu_run1/simu'
    mixer_simu = '/home/oleksandr/github/mixer/scripts/mixer_simu.py'
    tmp_prefix = '/scratch/simu_run1'

dry_run=False

num_submit = int(sys.argv[1])

basename_pattern = '_h2={h2}_rg={rg}_pi1u={pi1u}_pi2u={pi2u}_pifrac={pifrac}_rep={rep}'
cmd_pattern = 'bash && ~/miniconda3/bin/python3 {mixer_simu} --h2 {h2} --rg {rg} --pi1u {pi1u} --pi2u {pi2u} --pifrac {pifrac} --rep {rep} --tmp-prefix {tmp_basename} --out-prefix {out_basename}'

all_lists = []

pi1u_vec = [3e-4, 3e-3] # [3e-4, 3e-3, 3e-3]
pi2u_vec = [3e-4, 3e-3] # [3e-4, 3e-3, 3e-4]
pifrac_vec = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

for rep in range(1, 11):
    for h2 in ['0.4', '0.7', '0.1']: # ['0.4', '0.1', '0.7']:
        for rg in ['0.0', '0.5']:
            for pi1u, pi2u in zip(pi1u_vec, pi2u_vec):
                for pifrac in pifrac_vec:
                    all_lists.append([('rep', rep), ('h2', h2), ('rg', rg), ('pi1u', pi1u), ('pi2u', pi2u), ('pifrac', pifrac)])
del rep; del h2; del rg; del pi1u; del pi2u; del pifrac;

print(len(all_lists))

skipped = 0; duplicated = 0;
processed = set()  # exclude duplicates

cmd_list = []
touch_list = []
for val in all_lists:
    out_basename = out_prefix + basename_pattern.format(**locals(), **dict(val))
    if os.path.exists(out_basename + '.touch'): skipped+=1; continue

    tmp_basename = tmp_prefix + basename_pattern.format(**locals(), **dict(val))
    if out_basename in processed:
        duplicated+=1
        continue

    processed.add(out_basename)
    cmd = cmd_pattern.format(**locals(), **dict(val))

    analyses = ['simu', 'plink', 'ugmg', 'bgmg']
    
    cmd += ' --analysis ' + ' '.join(analyses)
    cmd_list.append(cmd)
    touch_list.append(out_basename+'.touch')

template='''
# Execute job in the queue "std.q" unless you have special requirements.
##$ -q std.q
##$ -q all_24.q

# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd

# Redirect output stream to this file.
##$ -o output.dat

# Redirect error stream to this file.
##$ -e error.dat

# Send status information to this email address.
##$ -M oleksandr.frei@gmail.com

# Send an e-mail when the job is done.
##$ -m e

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

bash
{}
'''

for cmd, touch in zip(cmd_list, touch_list):
    with open('run_script.sh', 'w') as f: f.write(template.format(cmd))
    if num_submit > 0:
        if not dry_run: os.system('qsub run_script.sh')
        if not dry_run: os.system('touch {}'.format(touch))
        print('qsub {}'.format(cmd))
    num_submit -= 1
    if num_submit <= 0: break

print('skipped: {}, duplicated: {}'.format(skipped, duplicated))

