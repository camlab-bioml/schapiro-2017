
import tifffile as tif
from scipy.io import loadmat
import os
import numpy as np
import pandas as pd

def parse_one_cell(image, mask, cell, channel, sample):
    is_cell = mask == cell
    x = np.where(is_cell)[0].mean()
    y = np.where(is_cell)[1].mean()
    size = is_cell.sum()
    expression = image[is_cell].mean()
    
    df = pd.DataFrame({
    'sample': [sample],
    'cell': [sample + "_" + str(cell)],
    'x': [x],
    'y': [y],
    'size': [size],
    'channel': [channel],
    'expression': [expression]
    })
    df = df.set_index("cell")
    return df
    
def parse_sample_single_channel(mask, image, channel, sample):
    cell_ids = np.unique(np.reshape(mask, -1))
    cell_ids = cell_ids[cell_ids != 0]
    pds = [parse_one_cell(image, mask, cell, channel, sample) for cell in cell_ids]
    return pd.concat(pds, axis=0)

def read_sample(mask, tiff_file, sample):
    channel = tiff_file.split("/")[-1].replace(".tiff", "")
    image = tif.imread(tiff_file)
    return parse_sample_single_channel(mask, image, channel, sample)

samples = pd.read_csv("samples.txt", header=None)
samples = list(samples[samples.columns[0]])

channels = pd.read_csv("channels.txt", header=None)
channels = list(channels[channels.columns[0]])




rule all:
    input:
        expand("single-cell-data/{s}.csv", s=samples),
        expand("single-cell-data/alt-masks/{alt_s}_{user}.csv", 
        alt_s = ['Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35'], user = ['Catena', 'Jackson', 'Schulz'])

rule to_csv:
    input:
        expand("raw-data/{{s}}/{{s}}/{channel}.tiff",
        channel=channels)
    output:
        "single-cell-data/{s}.csv"
    run:
        sample = wildcards.s
        sample_short = sample.split("_")[0]
        mask_file = f"raw-data/{sample}/CellProfiler/{sample_short}_Cells.mat"
        mask = loadmat(mask_file)
        mask = mask['Image']
        files = input


        df_list = [read_sample(mask, f, sample) for f in files]
        df = pd.concat(df_list, axis=0)

        df.to_csv(output[0])


rule to_csv_alt_masks:
    input:
        expand("raw-data/{{alt_s}}/{{alt_s}}/{channel}.tiff", channel=channels)
    output:
        "single-cell-data/alt-masks/{alt_s}_{user}.csv"
    run:
        sample = wildcards.alt_s
        user = wildcards.user
        mask_file = f"raw-data/alt_masks/{user}/{sample}_{user}_mask.tiff"
        mask = tif.imread(mask_file)
        files = input

        df_list = [read_sample(mask, f, sample) for f in files]
        df = pd.concat(df_list, axis=0)

        df.to_csv(output[0])

