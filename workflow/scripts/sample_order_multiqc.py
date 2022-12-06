#!/bin/python3

import sys #behovs?
import csv

samples = []
header = False

with open(snakemake.input.sample_sheet, 'r') as samplesheet:
    for lline in samplesheet:
        line = lline.strip()
        if header:
            samples.append(line.split(",")[1])
        if line == "Sample_ID,Sample_Name,Description,index,I7_Index_ID,index2,I5_Index_ID,Sample_Project":
            header = True


with open(snakemake.output[0], "w+") as target:
    tsv_writer = csv.writer(tsv, delimiter='\t')
    tsv_writer.writerow("Sample Name", "Sample Order")
    i = 1
    for sample in samples:
        tsv_writer.writerow(sample, 'sample'+str(i))
        i += 1
