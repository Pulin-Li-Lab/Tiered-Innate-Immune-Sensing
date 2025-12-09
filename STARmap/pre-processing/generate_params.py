# This script is executed by the qsub submission for steps 1, 2, and 3 of registration/spot-finding
# global registration and splitting of subtiles
# local registration and spot-finding
# stitching of subtiles

import os
from glob import glob
import pandas as pd
import argparse
import numpy as np

def generate_params(args):
    # Parameters to get: position number, sqrt_pieces, mode, subtile (1 for global/stitch), input_dim, sample
    sample = args.sample 
    mode = args.m
    if mode == 'global' or mode == 'local':
        mode = mode + '_registration'

    if args.tilelist:
        tilelist = args.tilelist
    else:
        tilelist = None
    
    if args.rounds:
        n_rounds = args.rounds
    else:
        n_rounds = 9
    
    if args.xy:
        xy = args.xy
    else:
        xy = 2048
    
    #samples = [f'{tissue}_STAR', f'{tissue}_RIBO'] 
    #tissue_dir = f'/stanley/WangLab/atlas/*{tissue}')[0]
    proj_dir = f'/stanley/WangLab/connie/{sample}'
    if not os.path.exists(f'{proj_dir}/00.sh/01.registration'):
        os.makedirs(f'{proj_dir}/00.sh/01.registration')

    with open(f'{proj_dir}/00.sh/01.registration/params_{mode.split("_")[0]}', 'w') as params_file:
        # Use sample name to get sample metadata
        meta = pd.read_csv('/stanley/WangLab/connie/sample_metadata.csv', index_col=0)
        z = meta.loc[sample, 'num_z']
        start = meta.loc[sample, 'start_tile']
        end = meta.loc[sample, 'end_tile']
        
        sqrt_pieces = 4
        input_dim = f'[{xy},{xy},{z},4,{n_rounds}]'
        subtile = 1

        if tilelist:
            positions = np.loadtxt(tilelist, dtype=int).flatten()
        else:
            positions = [str(pos).zfill(3) for pos in np.arange(start, end+1)]
           
        # Write into file
        for pos_num in positions:
            position = f'Position{pos_num}'
            if mode == 'local_registration':
                for subtile in range(1,sqrt_pieces**2 + 1):
                    params_file.write(f"{position};{sqrt_pieces};{mode};{subtile};{input_dim};{sample}\n")
            else:
                params_file.write(f"{position};{sqrt_pieces};{mode};{subtile};{input_dim};{sample}\n")
             

parser = argparse.ArgumentParser()
parser.add_argument("sample", type=str, help="sample name, ex. Heart")
parser.add_argument("-m", type=str, help="mode: global, local, or stitching")
parser.add_argument("-t", "--tilelist", type=str, help="optional path to list of Position numbers, for re-running on a subset of positions")
parser.add_argument("-r", "--rounds", type=int, help="optional number of rounds (defaults to 9)")
parser.add_argument("-xy", "--xy", type=int, help="optional x and y dimensions of tile (defaults to 2048)")
args = parser.parse_args()

generate_params(args)


    
