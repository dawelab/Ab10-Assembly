#!/usr/bin/env python3
import pandas as pd 
import numpy as np
import sys
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
import os
import subprocess
#from read_bionano import read_xmap,read_key


def read_xmap(xmap_file):
    xmap = pd.read_csv(xmap_file, delimiter="\t",comment="#",header=None).loc[:,[0,1,3,4,2,5,6,7,10,11,13]]
    xmap.columns = ['CMapId','Query_ID', 'Q_Start', 'Q_End','Ref_ID1', 'R_Start', 'R_End','Orientation','Q_Len', 'R_Len', 'Alignment_String']
    return xmap.reset_index(drop=True)
def read_key(key_file):
    key = pd.read_csv(key_file,delimiter="\t",comment="#",header=0)
    return key.reset_index(drop=True)

'''extract alignment block from xmap format'''
def extract_alignment(total):
    alignments = []    
    for x in total['Ref_ID1'].unique().tolist():
        refs = total[total['Ref_ID1']==x]
        #print(x)
        if len(refs['Query_ID'].unique().tolist()) >1: 
            for y in refs['Query_ID'].unique().tolist():
                querys = refs[refs['Query_ID']==y]
                orientation = max(set(querys['Orientation'].tolist()), key = querys['Orientation'].tolist().count)[0]
                r_start = min(querys['R_Start'].tolist())
                r_end = max(querys['R_End'].tolist())
                q_start = min(querys['Q_Start'].tolist() + querys['Q_End'].tolist())
                q_end = max(querys['Q_Start'].tolist() + querys['Q_End'].tolist())
                r_len = r_end - r_start
                q_len = abs(q_start - q_end)
                alignments.append([x,r_start,r_end,y,q_start,q_end,r_len,q_len,orientation])

    alignments = pd.DataFrame(alignments)    
    alignments.columns=['Ref_ID', 'R_Start','R_End','Query_ID', 'Q_Start','Q_End','R_Length','Q_Length','Orientation']
    return alignments
#alignments[alignments['R_Start']>=alignments['R_End']]
    
'''remove embedded contigs'''
def extract_embedded(alignments):
    embedded = pd.DataFrame()    
    for x in alignments['Ref_ID'].unique().tolist():
        refs = alignments[alignments['Ref_ID']==x].sort_values(by=['R_Start'])
        for i in refs.index:
            for j in refs.index:
                if refs.loc[i,'R_Start'] <= refs.loc[j,'R_Start'] and refs.loc[i,'R_End'] >= refs.loc[j,'R_End'] and i!=j:
                    embedded = embedded.append(refs.loc[j,:])
    return embedded       
 
'''identify overlaps between contigs, based on alignment to bionano maps'''
def identify_overlaps(alignments):
    overlaps = [] 
    for x in alignments['Ref_ID'].unique().tolist():
        refs = alignments[alignments['Ref_ID']==x].sort_values(by=['R_Start'])
        if len(refs.index) >1:
            for i in range(len(refs.index)-1):
                if refs.loc[refs.index[i+1],'R_Start'] < refs.loc[refs.index[i],'R_End']:
                    overlaps.append(refs.loc[refs.index[i],:].tolist() + refs.loc[refs.index[i+1],:].tolist())
    #$16 is overlapped sequence length
    overlaps = pd.DataFrame(overlaps)    
    overlaps.columns=['Ref1_ID', 'R1_Start','R1_End','Query1_ID', 'Q1_Start','Q1_End','R1_Length','Q1_Length','Orientation1','Ref2_ID', 'R2_Start','R2_End','Query2_ID', 'Q2_Start','Q2_End','R2_Length','Q2_Length','Orientation2']
    overlaps['Overlap_Len'] = (overlaps['R2_Start'] - overlaps['R1_End']).abs()
    #overlaps['Predicted_Merge_Len'] = (overlaps['R1_Length'] + overlaps['R2_Length'] - overlaps['Overlap_Len']).abs()
    return overlaps.sort_values('Ref1_ID').drop_duplicates(subset=['Query1_ID', 'Query2_ID'], keep='first').sort_index().reset_index(drop=True)

#map sequence ID with key
def map_key(overlaps,keys):
    for i in overlaps.index:
        r1key_id = overlaps.loc[i,'Query1_ID']
        r2key_id = overlaps.loc[i,'Query2_ID']
        if overlaps.loc[i,'R1_Length'] > overlaps.loc[i,'R2_Length']:
            overlaps.loc[i,'Ref_Len'] = overlaps.loc[i,'R1_Length']
            overlaps.loc[i,'Ref_Sequence_ID'] = keys[keys['CompntId']==r1key_id]['CompntName'].values[0]
            overlaps.loc[i,'Query_Len'] = overlaps.loc[i,'R2_Length']
            overlaps.loc[i,'Query_Sequence_ID'] = keys[keys['CompntId']==r2key_id]['CompntName'].values[0]
        else:
            overlaps.loc[i,'Ref_Len'] = overlaps.loc[i,'R2_Length']
            overlaps.loc[i,'Ref_Sequence_ID'] = keys[keys['CompntId']==r2key_id]['CompntName'].values[0]
            overlaps.loc[i,'Query_Len'] = overlaps.loc[i,'R1_Length']
            overlaps.loc[i,'Query_Sequence_ID'] = keys[keys['CompntId']==r1key_id]['CompntName'].values[0]
        if overlaps.loc[i,'Orientation1'] == overlaps.loc[i,'Orientation2']:
            overlaps.loc[i,'Orientation_final'] = "+"
        if overlaps.loc[i,'Orientation1'] != overlaps.loc[i,'Orientation2']:
            overlaps.loc[i,'Orientation_final'] = "-"
            
    return overlaps

#function to group
def group(overlaps_key):
    cluster = 0
    overlaps_key.loc[0,'cluster'] = 0
    for i in range(1,overlaps_key.shape[0]): 
        if overlaps_key.loc[i,'Ref1_ID'] == overlaps_key.loc[i-1,'Ref1_ID'] and overlaps_key.loc[i,'R1_Start'] <= overlaps_key.loc[i-1,'R1_End']:
            overlaps_key.loc[i,'cluster'] = cluster
            
        else:
            cluster += 1
            overlaps_key.loc[i,'cluster'] = cluster 
    return overlaps_key

#
'''determine ref and query in each overlap cluster'''
def ref_query(overlaps_key):
    for x in overlaps_key['cluster'].unique():
        overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==x]
        if overlaps_key_cluster.shape[0] == 1:
            i = overlaps_key_cluster.index[0]
            overlaps_key.loc[i,'Output_Sequence_ID'] =  'merge' + str(i) +'_map' + str(overlaps_key.loc[i,'Ref1_ID'])
            overlaps_key.loc[i,'Output_Predicted_Len'] = abs(overlaps_key.loc[i,'Ref_Len'] + overlaps_key.loc[i,'Query_Len'] - overlaps_key.loc[i,'Overlap_Len'])
            overlaps_key.loc[i,'Query_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Query_Sequence_ID']
            overlaps_key.loc[i,'Query_Len_final'] = overlaps_key_cluster.loc[i,'Query_Len']        
            overlaps_key.loc[i,'Ref_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Ref_Sequence_ID']
            overlaps_key.loc[i,'Ref_Len_final'] = overlaps_key_cluster.loc[i,'Ref_Len']
    
        #overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==4]
        if overlaps_key_cluster.shape[0] > 1:
            ref_query_list = []
    #        index0 = overlaps_key_cluster[overlaps_key_cluster[['Ref_Len']].values == overlaps_key_cluster[['Ref_Len']].values.max()].index.values[0] 
            i = overlaps_key_cluster.index[0]
            overlaps_key.loc[i,'Output_Sequence_ID'] =  'merge' + str(i) +'_map' + str(overlaps_key.loc[i,'Ref1_ID'])
            overlaps_key.loc[i,'Output_Predicted_Len'] = abs(overlaps_key.loc[i,'Ref_Len'] + overlaps_key.loc[i,'Query_Len'] - overlaps_key.loc[i,'Overlap_Len'])
            ref_query_list.extend([overlaps_key.loc[i,'Ref_Sequence_ID'],overlaps_key.loc[i,'Query_Sequence_ID']])
            overlaps_key.loc[i,'Query_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Query_Sequence_ID']
            overlaps_key.loc[i,'Query_Len_final'] = overlaps_key_cluster.loc[i,'Query_Len']        
            overlaps_key.loc[i,'Ref_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Ref_Sequence_ID']
            overlaps_key.loc[i,'Ref_Len_final'] = overlaps_key_cluster.loc[i,'Ref_Len']
            overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==x]
            for i in overlaps_key_cluster.index[1:]:
                overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==x]
                overlaps_key.loc[i,'Output_Sequence_ID'] = overlaps_key.loc[i-1,'Output_Sequence_ID'].split('_map')[0] +'_merge'+ str(i+1) +'_map' + str(overlaps_key.loc[i,'Ref1_ID'])
    #            overlaps_key.loc[i,'Output_Predicted_Len'] = abs(overlaps_key.loc[i,'Ref_Len'] + overlaps_key.loc[i,'Query_Len_final'] - overlaps_key.loc[i,'Overlap_Len'])
                if overlaps_key_cluster.loc[i,'Ref_Len'] > overlaps_key_cluster.loc[i-1,'Output_Predicted_Len']:
                    #print(i)
                    overlaps_key.loc[i,'Query_Sequence_ID_final'] = overlaps_key_cluster.loc[i-1,'Output_Sequence_ID']
                    overlaps_key.loc[i,'Query_Len_final'] = overlaps_key_cluster.loc[i-1,'Output_Predicted_Len']        
                    overlaps_key.loc[i,'Ref_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Ref_Sequence_ID']
                    overlaps_key.loc[i,'Ref_Len_final'] = overlaps_key_cluster.loc[i,'Ref_Len']
                    ref_query_list.append(overlaps_key.loc[i,'Ref_Sequence_ID'])
                else:
                    overlaps_key.loc[i,'Ref_Sequence_ID_final'] = overlaps_key_cluster.loc[i-1,'Output_Sequence_ID']
                    overlaps_key.loc[i,'Ref_Len_final'] = overlaps_key_cluster.loc[i-1,'Output_Predicted_Len']
                    query_ID =[q for q in overlaps_key.loc[i,['Ref_Sequence_ID','Query_Sequence_ID']].tolist() if not (q in ref_query_list)][0]
                    
                    if overlaps_key.loc[i,'Ref_Sequence_ID'] == query_ID:
                        #print(i)
                        #i=6
                        #print(overlaps_key.loc[i,'Ref_Sequence_ID'])
                        overlaps_key.loc[i,'Query_Sequence_ID_final'] = overlaps_key.loc[i,'Ref_Sequence_ID'] 
                        overlaps_key.loc[i,'Query_Len_final'] = overlaps_key.loc[i,'Ref_Len'] 
                        ref_query_list.append(overlaps_key.loc[i,'Ref_Sequence_ID'])
                    else:
                        overlaps_key.loc[i,'Query_Sequence_ID_final'] = overlaps_key_cluster.loc[i,'Query_Sequence_ID']
                        overlaps_key.loc[i,'Query_Len_final'] = overlaps_key_cluster.loc[i,'Query_Len']
                        ref_query_list.append(overlaps_key.loc[i,'Query_Sequence_ID'])
                overlaps_key.loc[i,'Output_Predicted_Len'] = abs(overlaps_key.loc[i,'Ref_Len_final'] + overlaps_key.loc[i,'Query_Len_final'] - overlaps_key.loc[i,'Overlap_Len'])
    return overlaps_key
#overlaps_key.loc[:,'Output_Sequence_ID']
#overlaps_key_cluster.loc[3,:]

'''output shell script for each overlap merging group, and execute shell'''
def output(overlaps_key,out_dir,genomepath):
    output_src_directory=os.path.join(out_dir, "src")
    os.makedirs(output_src_directory,exist_ok=True)
    modified_seq = []    
    output_seq = []
    for x in overlaps_key['cluster'].unique():
        overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==x]
        output_file='merge_overlap'+str(int(x))+'_map_'+str(overlaps_key[overlaps_key['cluster']==x]['Ref1_ID'].unique()[0])+'.sh'
        output_sh = open(os.path.join(output_src_directory, output_file), "w")
        output_sh.write('#!/bin/sh\n')
        overlap_sequence_ids = 'merge_overlap'+str(int(x))+'_map_'+str(overlaps_key[overlaps_key['cluster']==x]['Ref1_ID'].unique()[0])
        output_sh.write('mkdir -p ' + os.path.join(out_dir,overlap_sequence_ids) +'\n')
        output_sh.write('cd '+os.path.join(out_dir,overlap_sequence_ids) +'\n')   
        ref_query_all=overlaps_key_cluster[['Ref_Sequence_ID_final','Query_Sequence_ID_final']].values.tolist()
        ref_query=[item for sub in ref_query_all for item in sub if not 'merge' in item]            
        for item in ref_query:
            output_sh.write('samtools faidx ' + genomepath + ' ' + item + ' > ' + item + '.fasta\n')
        #print(os.path.join(output_src_directory, output_file))
        for i in overlaps_key_cluster.index:
            output_sh.write('cat ' + overlaps_key.loc[i,'Ref_Sequence_ID_final'] + '.fasta ' + overlaps_key.loc[i,'Query_Sequence_ID_final'] + '.fasta > ' + overlaps_key.loc[i,'Ref_Sequence_ID_final'] + '_' + overlaps_key.loc[i,'Query_Sequence_ID_final'] + '.fasta\n')
            output_sh.write('minimap2 -t22 -k28 -w28 -A1 -B9 -O16,41 -E2,1 -z200 -g100000 -r100000 --max-chain-skip 100 ' + overlaps_key.loc[i,'Ref_Sequence_ID_final'] + '.fasta ' + overlaps_key.loc[i,'Query_Sequence_ID_final'] + '.fasta > ' + overlaps_key.loc[i,'Output_Sequence_ID'] + '.paf\n')
            output_sh.write('cat ' + overlaps_key.loc[i,'Output_Sequence_ID'] + ".paf | awk '{if($5==" + '\"' + str(overlaps_key.loc[i,'Orientation_final']) + '\"' + "&&$11<=$2){print$0}}' > " + overlaps_key.loc[i,'Output_Sequence_ID'] + '.selected.paf\n')            
            output_sh.write('miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o' + str(int(overlaps_key.loc[i,'Overlap_Len'])) + ' ' + overlaps_key.loc[i,'Output_Sequence_ID'] + '.selected.paf > ' + overlaps_key.loc[i,'Output_Sequence_ID']  + '_noseq.gfa\n')
            output_sh.write('miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o' + str(int(overlaps_key.loc[i,'Overlap_Len'])) + ' ' + overlaps_key.loc[i,'Output_Sequence_ID'] + '.selected.paf -f ' + overlaps_key.loc[i,'Ref_Sequence_ID_final'] + '_' + overlaps_key.loc[i,'Query_Sequence_ID_final'] + '.fasta > ' + overlaps_key.loc[i,'Output_Sequence_ID'] + '.gfa\n')
            output_sh.write("awk '" + str('/^S/{print ">"$2"\\n"$3}') + "' " + overlaps_key.loc[i,'Output_Sequence_ID'] + '.gfa | fold > '+ overlaps_key.loc[i,'Output_Sequence_ID'] + '.fasta\n')
            output_sh.write("sed '" + str('s/>.*/&_') + overlaps_key.loc[i,'Output_Sequence_ID'] + "/g' " + overlaps_key.loc[i,'Output_Sequence_ID'] + '.fasta  > ' + overlaps_key.loc[i,'Output_Sequence_ID'] + '.renamed.fasta\n')  
        output_sh.close()
        os.chmod(os.path.join(output_src_directory, output_file), 0o775)
        subprocess.call([os.path.join(output_src_directory, output_file)])
        '''count the number of nucleotides in the output file'''
        for i in overlaps_key_cluster.index:    
            merged_file = os.path.join(overlap_sequence_ids, overlaps_key.loc[i,'Output_Sequence_ID']+ '.fasta')
            total=0
            with open(os.path.join(out_dir, merged_file), "r") as f:                
                for line in f:
                    if(line[0] == ">") :
                        continue
                    total += len(line)-1
            overlaps_key.loc[i,'Output_Len'] = total
            if overlaps_key.loc[i,'Output_Len']/overlaps_key.loc[i,'Output_Predicted_Len'] <= 1.5 and overlaps_key.loc[i,'Output_Len']/overlaps_key.loc[i,'Output_Predicted_Len'] >=0.5:
                modified_seq.extend([overlaps_key.loc[i,'Query_Sequence_ID_final'],overlaps_key.loc[i,'Ref_Sequence_ID_final']])

    for x in overlaps_key['cluster'].unique():
        overlaps_key_cluster = overlaps_key[overlaps_key['cluster']==x]    
        if overlaps_key_cluster.tail(1)['Output_Len'].values[0]/overlaps_key_cluster.tail(1)['Output_Predicted_Len'].values[0] <= 1.5 and overlaps_key_cluster.tail(1)['Output_Len'].values[0]/overlaps_key_cluster.tail(1)['Output_Predicted_Len'].values[0] >=0.5:
            output_seq.append(overlaps_key_cluster.tail(1)['Output_Sequence_ID'].values[0])
    modified_seq_input=[x for x in modified_seq if not 'merge' in x]
    overlaps_key.loc[:,['Ref_Sequence_ID_final','Ref_Len_final','Query_Sequence_ID_final','Query_Len_final',"Orientation_final",'cluster','Output_Sequence_ID','Output_Predicted_Len','Output_Len']].to_csv(os.path.join(out_dir,'overlaps_summary.txt'),sep='\t',float_format='%.1f')
    with open(os.path.join(out_dir,'modified_seq.txt'),'w') as f:
        f.writelines( "%s\n" % item for item in modified_seq_input)
    with open(os.path.join(out_dir,'output_seq.txt'),'w') as f:
        for item in output_seq:
            item_path = os.path.join(out_dir,'*',item+'.renamed.fasta')
            f.writelines(item_path + ' ')
    #print(output_seq)
       
def main(xmap_file,key_file,out_dir):
 
    total = read_xmap(xmap_file)
    #total.head()
    alignments = extract_alignment(total)
    embedded = extract_embedded(alignments)
    alignments.drop(set(embedded.index), axis=0,inplace=True)
    overlaps=identify_overlaps(alignments)
    keys = read_key(key_file)
    overlaps_key=map_key(overlaps,keys)
    overlaps_key=group(overlaps_key)
    overlaps_key_ref_query = ref_query(overlaps_key)
    genomepath=key_file.split('_CTTAAG_0kb_0labels_key.txt')[0] + '.fasta'
    output(overlaps_key_ref_query,out_dir,genomepath)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Detect overlaps among contigs after alignment to optical maps or hybrid scaffolding')
    parser.add_argument('-x','--xmap_file', required=True,
                        help='If hybrid assembly is performed, the xmap file is in directory HybridScaffold_output/hybrid_scaffolds; If only alignment is performed, the xmap file is in directory Alignment_output/alignref')
    parser.add_argument('-k','--key_file', required=True,
                        help='Key file of input fasta generated by fa2cmap_multi_color.pl, named _CTTAAG_0kb_0labels_key.txt')
    parser.add_argument('-o','--out_dir',required=True,
                        help='Output directory')
    args = parser.parse_args()
    main(xmap_file=args.xmap_file, key_file=args.key_file,out_dir=args.out_dir)

