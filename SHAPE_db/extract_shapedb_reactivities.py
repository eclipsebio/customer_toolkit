"""
Author: Dayanara Lebron-Aldea 
Release Date: 11/9/2022

This code intakes the SHAPE_{db}-{condition}_reactivities_{strand}.bedgraph file and intersects it with 
a region of interest provided in a BED file. It outputs a *.shape, *.fa and *.map that can be used as
input to folding algorithms. IQR or min-max scaling of reactivity values can be performed for outputs 
when the scaling parameter is set.
"""

import subprocess as sp
import argparse
import pandas as pd
import itertools
import os
import numpy as np

def read_bed(bedfile):
    """ Read a 6 column-based bed graph"""
    bed = pd.read_table(bedfile,header=None,names=["chr","start","stop","gene","num","strand"])
    strand = bed["strand"][0]
    return bed,strand

def IQR_scaling(series):
    """Scale raw reactivities with IQR"""
    IQR_factor = 1.5*(np.percentile(series,75)-np.percentile(series,25))
    scaled = round(series/IQR_factor,4)
    return scaled

def Min_Max_scaling(series):
    """Scale raw reactivities with Min-Max scaling"""
    min_point = float(series.min())
    max_point = float(series.max())
    scaling = max_point-min_point
    scaled = round((series - min_point)/scaling,4)
    return scaled

def output_map(df_bg,dfbed,scaling,gene_name):
    """subset bedgraph for region of interest and scale"""
    
    exon_subsets = list() 
    df = pd.DataFrame()
    pos = list()
    
    #Make sure these are int and not float;
    dfbed.loc[:,'stop'] = dfbed.loc[:,'stop'].astype(int)
    dfbed.loc[:,'start'] = dfbed.loc[:,'start'].astype(int)
    
    dfbed.loc[:,'length'] = dfbed['stop'] - dfbed['start']
    total_length = dfbed['length'].sum()
    
    chrom = dfbed['chr'].unique()[0]
    
    df_chrom = df_bg[df_bg['chr'] == chrom]
    
    dfbed.loc[:,"length"] = dfbed["stop"]-dfbed["start"]
    
    #Better to do the merge after the scalings to input the -999
    for i,exon in dfbed.iterrows():
        start = int(exon['start'])
        stop = int(exon['stop'])
        pos.append(list(range(start,stop)))
        exon_subsets.append(df_chrom[df_chrom['start'].between(start,stop-1)])
                           
    df.loc[:,'start'] = list(itertools.chain(*pos))
    df.loc[:,'stop'] = df['start'] + 1
    
    #Reactivities for all exons, concatenated as rows
    region_map = pd.concat(exon_subsets,axis=0)
    
    if scaling == 'IQR':
        region_map.loc[:,'scaled'] = IQR_scaling(region_map['mean'].astype(float))
        print(f"Computing {scaling} scaled reactivities")                         
    elif scaling == 'Min-Max':
        region_map.loc[:,'scaled'] = Min_Max_scaling(region_map['mean'].astype(float))
        print(f"Computing {scaling} scaled reactivities") 
                           
    
    #Fill in the gaps with -999
    region_map_filled = df.merge(region_map,on=['start','stop'],how='outer').replace(np.nan,-999) 
  
    region_map_filled.loc[:,'start'] = region_map_filled.loc[:,'start'].astype(int)
    region_map_filled.loc[:,'stop'] = region_map_filled.loc[:,'stop'].astype(int) 
    positions = list(range(1,total_length+1))
    region_map_filled.loc[:,'1-based'] = positions
                           
    if region_map_filled.shape[0] != total_length:
        print(f'Selected dataframe doesnt equal region length of {total_length}; DF Shape is {region_map.shape[0]}')  
    
    region_map_filled.loc[:,'chr'] = [chrom]*total_length  
    
    if scaling:
        shape_cols= ['1-based','scaled']
        map_cols = ['chr','start','stop','mean','scaled']
        out_map = f'{gene_name}_{scaling}.map'
    else:
        shape_cols = ['1-based','mean']
        map_cols = ['chr','start','stop','mean']
        out_map = f'{gene_name}.map'
        
    df_shape = region_map_filled[shape_cols]                 
    region_map_filled[map_cols].to_csv(out_map,sep="\t",index = None)
    df_shape.to_csv(out_map.replace('map','shape'),header = None, sep="\t",index=None)
    
    forna_reactivity_gradient(df_shape,f'{gene_name}_forna_colors.txt')
    
    print('Finished')
                     
def read_bedgraph(bedgraph):
    """Read map file and output gene fasta"""
    df = pd.read_table(bedgraph,header = None,sep = ' ',names = ['chr','start','stop','mean'],skipfooter=1, engine = 'python')
    
    # Replace  -999 (positions w/o reactivities) to NAN so it doesnt affect calculations
    df.replace(-999,np.nan,inplace = True)
    
    return df

def intersect_fasta(reference,bed):
    gene = bed.replace('.6.bed','').replace('.bed','')
    fasta_output = gene +".fasta"
    tmp_fasta = fasta_output.replace('.fasta','.tmp')
    cmd = f'bedtools getfasta -fi {reference} -bed {bed} -fo {tmp_fasta}'
    sp.call(cmd , shell = True)
    
    ## Get intersected sequence to format accepted by algorithms
    out = open(fasta_output, 'w')
    out.write(f'>{gene}\n')
    with open(tmp_fasta) as t:
        for line in t.readlines():
            if not line.startswith('>'):
                out.write(line.strip().replace('T','U'))
    out.close()
    os.remove(tmp_fasta)
    
    print(f'Intersected Fasta is available here: {fasta_output}')
    
def forna_reactivity_gradient(df_shape,output):
    out = open(output,'w')
    df_forna = df_shape.copy()
    if 'scaled' in df_forna.columns:
        reactivity_column = 'scaled'
    else:
        reactivity_column = 'mean'
    
    df_forna.loc[:,"na"] = df_forna.loc[:,reactivity_column].replace(-999,np.nan)

    r_min = df_forna['na'].min(skipna=True).round(4)
    r_max = df_forna['na'].max(skipna=True).round(4)
    
    #Write forna_colors.txt
    out.write("range=aqua:red\n")
    out.write(f"domain = {r_min}:{r_max}\n")
    for element in list(df_forna[reactivity_column]):
        out.write(str(element) + "\n")
    out.close()
    


if __name__ == "__main__":

    help_string = "This script intersects the input bedgraph for a region of interest. The tool assume consecutive/non-overlapped exon locations in the bed file. Scaling options are: 'IQR and Min-Max'; Outputs: .shape and fasta files together with a forna_colors.txt"
        
    parser=argparse.ArgumentParser(description = help_string)
    parser.add_argument('--bed',help = 'BED file for region of interest', required = True)
    parser.add_argument('--bedgraph',help = 'Stranded bedgraph file',required=True)
    parser.add_argument('--fasta', help = 'Reference fasta to extract sequence from (optional)',default=None)
    parser.add_argument('--scaling', help='Scaling options: "IQR" or "Min-Max"; Otherwise outputs raw reactivities')
    args=parser.parse_args()
    
    strand_to_str = {'pos':'+','neg':'-'}
    target = os.path.basename(args.bed).split(".")[0] #grab basename
    
    print(f'Analyzing target: {target}')
    
    #check that the bedgraph ends in bedgraph
    basename = args.bedgraph.replace(".bedgraph","")
    strand = basename.split("_")[-1] #this will say pos and neg
                        
    #Check that strand information matches bed information
    dfbed,strandbed = read_bed(args.bed)
    
    if strandbed == strand:
        pass
    elif strandbed != strand_to_str[strand]:
        print("Bed file strand does not match bedgraph strand, please provide correct bedgraph")
        exit()
               
    if len(dfbed.columns) != 6:
        print("Please provide a 6 column bed")
        exit()
    
    if (args.scaling != None) and (args.scaling not in ['IQR','Min-Max']):
        print('Options for scaling are: IQR or Min-Max')
        exit()
    
    if len(dfbed['chr'].unique()) > 1:
        print('Bed File may only contain regions from 1 chromosome.')
        exit()
        
        
    #Check that bedfile has non-overlapping exon
    df_bedgraph = read_bedgraph(args.bedgraph)

    ## Subset the sequence in the bed region from a reference.
    if args.fasta:
        intersect_fasta(args.fasta,args.bed)
    
    output_map(df_bedgraph,dfbed,args.scaling,target)
