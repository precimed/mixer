#!/bin/bash

cfg_templ=optimize.e1.template.cfg

for i_hsq in $(seq 0 2); do
    for i_coding in $(seq 0 2); do
        for i_noncoding in $(seq 0 2); do
            for i_s2coding in $(seq 0 2); do
                for i_repeat in $(seq 10); do
                    cfg_f=optimize.e1.hsq_${i_hsq}.coding_${i_coding}.noncoding_${i_noncoding}.s2coding_${i_s2coding}.${i_repeat}.cfg
                    cp ${cfg_templ} ${cfg_f}
                    sed -i "s/i_hsq = 0/i_hsq = ${i_hsq}/" ${cfg_f}
                    sed -i "s/i_coding = 0/i_coding = ${i_coding}/" ${cfg_f}
                    sed -i "s/i_noncoding = 0/i_noncoding = ${i_noncoding}/" ${cfg_f}
                    sed -i "s/i_s2coding = 0/i_s2coding = ${i_s2coding}/" ${cfg_f}
                    sed -i "s/i_repeat = 1/i_repeat = ${i_repeat}/" ${cfg_f}
                    echo ${cfg_f} created
                done
            done
        done
    done
done
