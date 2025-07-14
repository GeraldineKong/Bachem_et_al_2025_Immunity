#!/usr/bin/env python

import pandas as pd 
import argparse

def merge_files(checkm2, comp, cont, gunc_res):
    #file = pd.read_csv('/home/kongg1/shotgun_data_2023/mouse_single_assembly/bins_vamb_single_drep_checkm/quality_report.tsv', sep = "\t")
    file = pd.read_csv(checkm2,sep = "\t")
    file = file[(file.Completeness >= comp) & (file.Contamination <= cont)]

    #gunc_file = pd.read_csv('/home/kongg1/shotgun_data_2023/mouse_single_assembly/bin_vamb_single_drep_gunc/GUNC.progenomes_2.1.maxCSS_level.tsv', sep = "\t")
    gunc_file = pd.read_csv(gunc_res, sep = "\t")
    gunc_file = gunc_file[gunc_file['pass.GUNC'] == True]

    file_final = file.merge(gunc_file, how = 'inner', left_on = 'Name', right_on = 'genome')
    return(file_final)
    

def main():
    parser = argparse.ArgumentParser(description = 'Find bins that pass checkm2 and GUNC QC')
    parser.add_argument('--checkm2', help = 'path/to/checkm2/res/quality_report.tsv')
    parser.add_argument('--comp', default = 90, type = int, help = 'Completion cutoff for checkm2. (Default: 90)')
    parser.add_argument('--cont', default = 5, type = int, help = 'Contamination cutoff for checkm2. (Default: 5)')
    parser.add_argument('--gunc', help = 'path/to/gunc/res/GUNC.progenomes_2.1.maxCSS_level.tsv')
    parser.add_argument('--output', help = 'path/to/output.txt')

    args = parser.parse_args()

    # Merge files
    file_final = merge_files(args.checkm2, args.comp, args.cont, args.gunc)
    
    # Save output
    file_final.to_csv(args.output, sep = "\t", index = False)

if __name__ == "__main__":
    main()