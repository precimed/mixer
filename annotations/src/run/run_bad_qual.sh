#!/bin/bash

bad_qual_f="bad_qual_070519.csv"

while read i_h2 i_p_c i_p_nc i_s2 i_repeat; do
    sbatch run_optimize.e1.sh ${i_h2} ${i_p_c} ${i_p_nc} ${i_s2} ${i_repeat};
done < ${bad_qual_f}

