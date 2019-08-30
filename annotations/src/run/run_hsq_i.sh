#!/bin/bash

i_hsq=2

for i_coding in $(seq 0 2); do
    for i_noncoding in $(seq 0 2); do
        for i_s2coding in $(seq 0 2); do
            for i_repeat in $(seq 10); do
                sbatch run_optimize.e1.sh ${i_hsq} ${i_coding} ${i_noncoding} ${i_s2coding} ${i_repeat}
            done
        done
    done
done
 
