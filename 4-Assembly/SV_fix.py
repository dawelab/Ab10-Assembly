
import pandas as pd 
import numpy as np
import sys
import re
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
import glob
from read_bionano import read_smap,read_xmap,read_cmap,read_key
'''
Files required:
    1) Smap file from PacBio
    2) Smap file from Nanopore
    4) Xmap from Nanopore
    6) Cmap from Nanopore
    7) Key file from PacBio
    8) Key file from Nanopore
'''

'''Identify SVs unique in PacBio assembly by comparing Reference coordinates in SV outputs for two assemblies, 
   Output Ref coordinates where SVs were found in PacBio'''
def SV_unique(smap_pbio_pav,smap_nanopore_pav):
    refIDs = smap_pbio_pav['Ref_ID1'].unique().tolist()
    pav_overlap = []
    for x in refIDs:
        smap_pbio_pav_ID = smap_pbio_pav[smap_pbio_pav['Ref_ID1'] == x]
        smap_nanopore_pav_ID = smap_nanopore_pav[smap_nanopore_pav['Ref_ID1'] == x]
        for i in smap_pbio_pav_ID.index:
            for j in smap_nanopore_pav_ID.index:
                if len(smap_nanopore_pav_ID.index) >0:
                    if smap_nanopore_pav_ID.loc[j,'R_Start'] >= smap_pbio_pav_ID.loc[i,'R_Start'] and smap_nanopore_pav_ID.loc[j,'R_Start'] <= smap_pbio_pav_ID.loc[i,'R_End']:
                        pav_overlap.append(i)
                    elif smap_nanopore_pav_ID.loc[j,'R_End'] >= smap_pbio_pav_ID.loc[i,'R_Start'] and smap_nanopore_pav_ID.loc[j,'R_End'] <= smap_pbio_pav_ID.loc[i,'R_End']:
                        pav_overlap.append(i)
                    elif smap_nanopore_pav_ID.loc[j,'R_Start'] <= smap_pbio_pav_ID.loc[i,'R_Start'] and smap_nanopore_pav_ID.loc[j,'R_End'] >= smap_pbio_pav_ID.loc[i,'R_End']:
                        pav_overlap.append(i)
    pav_ref_replace_Index = set(smap_pbio_pav.index) - set(pav_overlap)    
    return smap_pbio_pav.loc[set(pav_ref_replace_Index)].loc[:,['Ref_ID1','R_Start', 'R_End','Query_ID','Q_Start','Q_End','RefStartIdx','RefEndIdx','QryStartIdx', 'QryEndIdx','XmapID1','SV_Type','SV_Size']]


'''Use PacBio unique SV coordinates as input, identify corresponding coordinates and Index in Nanopore xmap file where assemblies are accurate '''

def find_idx(RefIdx_PacBio,xmap_nanopore):
    nanopore_Idx_all_start = []
    nanopore_Idx_all_end = []    
    for x in RefIdx_PacBio['Ref_ID1'].unique().tolist():
        RefIdx_PacBio_IDs = RefIdx_PacBio[RefIdx_PacBio['Ref_ID1']==x]
        xmap_nanopore_IDs = xmap_nanopore[xmap_nanopore['Ref_ID1']==x]
        #xmap_nanopore_IDs.head()
        #RefIdx_PacBio_IDs.head()
        for num in RefIdx_PacBio_IDs.index: 
            startsIdx = RefIdx_PacBio_IDs.loc[num,'RefStartIdx']
            startsPos = RefIdx_PacBio_IDs.loc[num,'R_Start']
            endsIdx = RefIdx_PacBio_IDs.loc[num,'RefEndIdx']
            endsPos = RefIdx_PacBio_IDs.loc[num,'R_End']
            pacbio_ID = RefIdx_PacBio_IDs.loc[num,'Query_ID']
            pacbio_startsIdx = RefIdx_PacBio_IDs.loc[num,'QryStartIdx']
            pacbio_startsPos = RefIdx_PacBio_IDs.loc[num,'Q_Start']
            pacbio_endsIdx = RefIdx_PacBio_IDs.loc[num,'QryEndIdx']
            pacbio_endsPos = RefIdx_PacBio_IDs.loc[num,'Q_End']
            pacbio_svtype = RefIdx_PacBio_IDs.loc[num,'SV_Type']
            pacbio_svsize = RefIdx_PacBio_IDs.loc[num,'SV_Size']
        #    starts_var = [startsIdx-1,startsIdx,startsIdx+1]
            for i in xmap_nanopore_IDs.index:
                nanopore_IDs_Query_ID = xmap_nanopore_IDs.loc[i,'Query_ID']
                if startsPos < xmap_nanopore_IDs.loc[i,'R_End'] and startsPos > xmap_nanopore_IDs.loc[i,'R_Start']:
                    #print(i)
                    alignment_string = re.findall("[^()]+", xmap_nanopore_IDs.loc[i,'Alignment_String'])
                    alignment_string_reformat = pd.DataFrame([j.split(",") for j in alignment_string])
                    alignment_string_reformat.columns = ['RefIdx','QueryIdx']
                    for alignment_index in alignment_string_reformat.index:
                        if alignment_string_reformat.loc[alignment_index,'RefIdx'] == str(startsIdx):
                            #print(i)
                            nanopore_Idx_Start = int(''.join(alignment_string_reformat[alignment_string_reformat['RefIdx'] == str(startsIdx)]['QueryIdx'].tolist()))
                            nanopore_Idx_all_start.append([num,x,startsIdx,startsPos,pacbio_ID,pacbio_startsIdx,pacbio_startsPos,nanopore_IDs_Query_ID,nanopore_Idx_Start,pacbio_svtype,pacbio_svsize])
    
                if endsPos < xmap_nanopore_IDs.loc[i,'R_End'] and endsPos > xmap_nanopore_IDs.loc[i,'R_Start']:
                    #print(i)
                    alignment_string = re.findall("[^()]+", xmap_nanopore_IDs.loc[i,'Alignment_String'])
                    alignment_string_reformat = pd.DataFrame([j.split(",") for j in alignment_string])
                    alignment_string_reformat.columns = ['RefIdx','QueryIdx']
                    for alignment_index in alignment_string_reformat.index:
                        if alignment_string_reformat.loc[alignment_index,'RefIdx'] == str(endsIdx):
                            #print(i)
                            nanopore_Idx_end = int(''.join(alignment_string_reformat[alignment_string_reformat['RefIdx'] == str(endsIdx)]['QueryIdx'].tolist()))
                            nanopore_Idx_all_end.append([num,x,endsIdx,endsPos,pacbio_ID,pacbio_endsIdx,pacbio_endsPos,nanopore_IDs_Query_ID,nanopore_Idx_end,pacbio_svtype,pacbio_svsize])
    nanopore_Idx_all_start = pd.DataFrame(nanopore_Idx_all_start)
    nanopore_Idx_all_start.columns=['SV_index','BioNano_ID','BioNano_startsIdx', 'BioNano_startsPos','PacBio_ID','PacBio_startsIdx', 'PacBio_startsPos','Nanopore_ID','Nanopore_startsIdx', 'SV_Type','SV_Size']
    nanopore_Idx_all_end = pd.DataFrame(nanopore_Idx_all_end)
    nanopore_Idx_all_end.columns=['SV_index','BioNano_ID','BioNano_endsIdx', 'BioNano_endsPos','PacBio_ID','PacBio_endsIdx', 'PacBio_endsPos','Nanopore_ID','Nanopore_endsIdx', 'SV_Type','SV_Size']
    return nanopore_Idx_all_start,nanopore_Idx_all_end

'''Find corresponding coordinates in cmap with Idx'''
def map_coordinates(nanopore_paired_Idx,cmap_nanopore):
    nanopore_paired_Idx.insert(20, "Nanopore_startsPos", 'NA') 
    nanopore_paired_Idx.insert(21, "Nanopore_endsPos", 'NA') 
    for i in nanopore_paired_Idx.index:
        nanopore_ID = nanopore_paired_Idx.loc[i,'Nanopore_ID']
        nanopore_startsIdx = nanopore_paired_Idx.loc[i,'Nanopore_startsIdx']
        nanopore_endsIdx = nanopore_paired_Idx.loc[i,'Nanopore_endsIdx']
        nanopore_paired_Idx.loc[i,'Nanopore_startsPos'] = cmap_nanopore[(cmap_nanopore['CMapId']==nanopore_ID)&(cmap_nanopore['SiteID']==nanopore_startsIdx)]['Position'].values[0]
        nanopore_paired_Idx.loc[i,'Nanopore_endsPos'] = cmap_nanopore[(cmap_nanopore['CMapId']==nanopore_ID)&(cmap_nanopore['SiteID']==nanopore_endsIdx)]['Position'].values[0]
    return nanopore_paired_Idx


'''Find corresponding sequence ID in key file'''
def map_sequence_IDs(nanopore_keys,pacbio_keys,pav_sum):
    pav_sum.insert(22, "PacBio_sequence_ID", 'NA') 
    pav_sum.insert(23, "Nanopore_sequence_ID", 'NA') 
    
    for i in pav_sum.index:
        pacbio_ID = pav_sum.loc[i,'PacBio_ID_x']
        nanpore_ID = pav_sum.loc[i,'Nanopore_ID']
        pav_sum.loc[i,'PacBio_sequence_ID'] = pacbio_keys[pacbio_keys['CompntId']==pacbio_ID]['CompntName'].values[0]
        pav_sum.loc[i,'Nanopore_sequence_ID'] = nanopore_keys[nanopore_keys['CompntId']==nanpore_ID]['CompntName'].values[0]
    return pav_sum

'''Extract sequence by ID from multiple fasta file'''
def extract_seq(seq_record, ID):
    for i in range(len(seq_record)):
        if seq_record[i].id == str(ID):
            return seq_record[i]

def fix_SV_seq(seq_2_modify,position_list,seq_record_nanopore):
    modified_seq = seq_2_modify.seq[:int(position_list.loc[0,'PacBio_startsPos'])]
#len(modified_seq)
    for i in position_list.index:
        #print(i)
        seq_supply = extract_seq(seq_record_nanopore, position_list.loc[i,'Nanopore_sequence_ID'])
        if position_list.loc[i,'Orientation'] == "+":    
            modified_seq += seq_supply.seq[int(position_list.loc[i,'Nanopore_startsPos']):int(position_list.loc[i,'Nanopore_endsPos'])]
        else:
            modified_seq += seq_supply.seq[int(position_list.loc[i,'Nanopore_startsPos']):int(position_list.loc[i,'Nanopore_endsPos'])].reverse_complement()
        if i < max(position_list.index):
            #print(i)
            modified_seq += seq_2_modify.seq[int(position_list.loc[i,'PacBio_endsPos']):int(position_list.loc[i+1,'PacBio_startsPos'])]
        else:
            #print(i)
            modified_seq += seq_2_modify[int(position_list.loc[i,'PacBio_endsPos']):]
    return modified_seq

'''Generate list of coordinates and sequences from both PacBio and Nanopore fasta files for merging'''
def rewrite_records(pav_sum_sequence,seq_record_pacbio,seq_record_nanopore):
    all_records = []
    # Extract sequence of IDs with fixable SV
    for x in pav_sum_sequence['PacBio_sequence_ID'].unique().tolist():
        #print(x)
        pav_sum_sequence_ID  = pav_sum_sequence[pav_sum_sequence['PacBio_sequence_ID']==x].sort_values(by=['PacBio_startsPos']).reset_index(drop=True)
        seq_record_pacbio_ID = extract_seq(seq_record_pacbio, x)
        revised_seq = fix_SV_seq(seq_record_pacbio_ID,pav_sum_sequence_ID,seq_record_nanopore)
        all_records.append(revised_seq)       
    seq_record_pacbio_allIDs = []
    # Extract sequence of Ids with no PAV
    for i in range(len(seq_record_pacbio)):
        seq_record_pacbio_allIDs.append(seq_record_pacbio[i].id)
    unmodified_IDs = set(seq_record_pacbio_allIDs) - set(str(i) for i in pav_sum_sequence['PacBio_sequence_ID'].unique().tolist())
    #len(unmodified_IDs)
    for ids in unmodified_IDs:
        all_records.append(extract_seq(seq_record_pacbio,ids))
    all_records.sort(key=lambda r:int(r.id))
    return all_records



def main(core_SV_dir,sup_SV_dir,core_cut_dir,sup_cut_dir,output_dir,SV_size_cutoff):
    
    pacbio_smap_file = glob.glob(core_SV_dir +'/*_0kb_0labels.smap')[0]
    nanopore_smap_file = glob.glob(sup_SV_dir +'/*_0kb_0labels.smap')[0]
    nanopore_xmap_file = glob.glob(sup_SV_dir +'/*_0kb_0labels.xmap')[0]
    nanopore_cmap_file = glob.glob(sup_SV_dir +'/*_0kb_0labels.cmap')[0]
    pacbio_key_file = glob.glob(core_cut_dir +'/*_0kb_0labels_key.txt')[0]
    pacbio_fasta_file = glob.glob(core_cut_dir +'/*_cut.fasta')[0] 
    nanopore_key_file = glob.glob(sup_cut_dir +'/*_0kb_0labels_key.txt')[0]    
    nanopore_fasta_file = glob.glob(sup_cut_dir +'/*_cut.fasta')[0] 
    output_fasta_file = output_dir + '/core_SV_fixed.fasta'
    output_svcoord_file = output_dir + '/SV_fix_coord_sum.txt'
    
    '''Read smap files -- SV files, xmap files -- alignment files, of PacBio and Nanopore assemblies'''
    smap_pbio_pav = read_smap(pacbio_smap_file)[0]
    smap_nanopore_pav = read_smap(nanopore_smap_file)[0]
    #smap_nanopore_pav[(smap_nanopore_pav['Ref_ID1'] == 6)&(smap_nanopore_pav['Query_ID']==64)]
    xmap_nanopore = read_xmap(nanopore_xmap_file)
    ''' Read cmap file from Nanopore, and key file of PacBio and Nanopore assemblies (used to identify coordinates with index)'''
    cmap_nanopore = read_cmap(nanopore_cmap_file)
    nanopore_keys = read_key(nanopore_key_file)
    pacbio_keys = read_key(pacbio_key_file)
    
    '''Identify indels present in PacBio assembly, where no SVs were found in Nanopore assembly (correct assembly due to long reads)'''
    RefIdx_PacBio = SV_unique(smap_pbio_pav,smap_nanopore_pav)
    
    '''Number and size of total indels identified in PacBio assembly'''
    print(str(smap_pbio_pav.shape[0]) + ' assembly errors (>=1kb) are present in the core assembly.')
    print('Among the assembly errors, ' + str(smap_pbio_pav[smap_pbio_pav['SV_Type'].str.contains("deletion")].shape[0]) + ' are deletion errors, with total size of ' + str(round(smap_pbio_pav[smap_pbio_pav['SV_Type'].str.contains("deletion")]['SV_Size'].sum()/1000000,2)) +' Mb; ' + \
    str(smap_pbio_pav[smap_pbio_pav['SV_Type'].str.contains("insertion")].shape[0]) + ' are insertion errors, with total size of  ' + str(round(smap_pbio_pav[smap_pbio_pav['SV_Type'].str.contains("insertion")]['SV_Size'].sum()/1000000,2)) +' Mb.') 
    '''Number and size of indels unique in PacBio assembly'''
    print(str(RefIdx_PacBio.shape[0]) + ' out of ' + str(smap_pbio_pav.shape[0]) + ' misassembly errors are specific in the core assembly, and could potentially be fixed by the supplementary assembly.')
    print(str(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("deletion")].shape[0]) + ' deletion errors of ' + str(round(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("deletion")]['SV_Size'].sum()/1000000,2)) +' Mb in total. ' + \
    str(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("insertion")].shape[0]) + ' insertion erros of ' + str(round(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("insertion")]['SV_Size'].sum()/1000000,2)) +' Mb in total.') 
    
    #RefIdx_PacBio.sort_values(by=['SV_Size']).tail(10)
    '''Use Bionano as an anchor, identify the nanopore assembly coordinates where PacBio assembly has indels'''
    
    nanopore_startIdx = find_idx(RefIdx_PacBio,xmap_nanopore)[0]
    nanopore_startIdx = nanopore_startIdx.sort_values(by=['Nanopore_ID']).drop_duplicates(subset=['SV_index'],keep='first')
    nanopore_endIdx = find_idx(RefIdx_PacBio,xmap_nanopore)[1]
    nanopore_endIdx = nanopore_endIdx.sort_values(by=['Nanopore_ID']).drop_duplicates(subset=['SV_index'],keep='first')
    
    nanopore_paired_Idx = pd.merge(nanopore_startIdx, nanopore_endIdx, on=['SV_index','Nanopore_ID'],how='inner')
    print('Among the ' + str(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("deletion")].shape[0]) + ' deletion errors, ' + str(nanopore_paired_Idx[nanopore_paired_Idx['SV_Type_y'].str.contains("deletion")].shape[0]) +' will be fixed by the supplementary assembly due to its high-confidence alignment with optical maps, of size ' +str(round(nanopore_paired_Idx[nanopore_paired_Idx['SV_Type_y'].str.contains("deletion")]['SV_Size_y'].sum()/1000000,2)) +' Mbp.')
    print('Among the ' + str(RefIdx_PacBio[RefIdx_PacBio['SV_Type'].str.contains("insertion")].shape[0]) + ' insertion errors, ' + str(nanopore_paired_Idx[nanopore_paired_Idx['SV_Type_y'].str.contains("insertion")].shape[0]) +' will be fixed by the supplementary assembly due to its high-confidence alignment with optical maps, of size ' +str(round(nanopore_paired_Idx[nanopore_paired_Idx['SV_Type_y'].str.contains("insertion")]['SV_Size_y'].sum()/1000000,2)) +' Mbp.')
    
    nanopore_paired_Idx = nanopore_paired_Idx[nanopore_paired_Idx['SV_Size_y']>=SV_size_cutoff]
    
    
    '''Identify coordinates and ID in fasta, where sequence replacement occur'''
    pav_sum = map_coordinates(nanopore_paired_Idx,cmap_nanopore)
    pav_sum = map_sequence_IDs(nanopore_keys,pacbio_keys,pav_sum)
    pav_sum_sequence = pav_sum.loc[:,['PacBio_sequence_ID','PacBio_startsPos','PacBio_endsPos','Nanopore_sequence_ID','Nanopore_startsPos','Nanopore_endsPos']]
    pav_sum_sequence['Orientation'] = ['+' if pav_sum_sequence.loc[i,'Nanopore_startsPos'] < pav_sum_sequence.loc[i,'Nanopore_endsPos'] else '-' for i in pav_sum_sequence.index]
    
    
    '''Fasta files from PacBio and Nanopore contigs, for sequence replacement'''
    seq_record_pacbio = list(SeqIO.parse(pacbio_fasta_file, "fasta"))
    seq_record_nanopore = list(SeqIO.parse(nanopore_fasta_file, "fasta"))
    pacbio_records = rewrite_records(pav_sum_sequence,seq_record_pacbio,seq_record_nanopore)
    #len(pacbio_records)
    SeqIO.write(pacbio_records, output_fasta_file, "fasta")
    pav_sum_sequence.rename(columns={"PacBio_sequence_ID": "core_sequence_ID", "PacBio_startsPos": "core_startsPos","PacBio_endsPos": "core_endsPos", "Nanopore_sequence_ID": "supplementary_sequence_ID","Nanopore_startsPos": "supplementary_startsPos","Nanopore_endsPos": "supplementary_endsPos","Orientation": "relative_orientation"})
    pav_sum_sequence.to_csv(output_svcoord_file,sep='\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Cut fasta files based on conflict junction coordinates on cmap level')
    parser.add_argument('-cs','--core_SV_dir', default='./core_assembly/step2_SV_fix/SVs_identified',
                        help='SV output director for the core genome')
    parser.add_argument('-ss','--sup_SV_dir', default='./sup_assembly/step2_SV_fix/SVs_identified',
                        help='SV output director for the supplementary genome')
    parser.add_argument('-cc','--core_cut_dir',default='./core_assembly/step1_cut_conflicts/conflict_NGS',
                        help='auto_cut_NGS_coord_translation.txt file in directory cut_conflicts, generated by cut_conflict.pl')
    parser.add_argument('-sc','--sup_cut_dir', default='./sup_assembly/step1_cut_conflicts/conflict_NGS',
                        help='auto_cut_NGS_coord_translation.txt file in directory cut_conflicts, generated by cut_conflict.pl')
    parser.add_argument('-o','--output_dir', default='./core_assembly/step2_SV_fix',
                        help='output fasta file with conflict resolved on the sequence level')
    parser.add_argument('-s','--SV_size_cutoff', default=1000,
                        help='Misassemblies in core assembly above this threshold will be fixed')

    args = parser.parse_args()
    main(core_SV_dir=args.core_SV_dir, sup_SV_dir=args.sup_SV_dir,core_cut_dir=args.core_cut_dir, sup_cut_dir=args.sup_cut_dir,output_dir=args.output_dir,SV_size_cutoff=args.SV_size_cutoff)





