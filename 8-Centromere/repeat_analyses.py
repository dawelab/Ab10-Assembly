#!/usr/bin/env python3

import pandas as pd 
import os
import numpy as np

'''Read blast output, convert coordinates as start < end, and mark forward/reverse orientation'''
def reorder(blast_file):
    blast = pd.read_csv(blast_file, delimiter="\t",header=None)
    blast_direct_orientation = blast[blast[8]<blast[9]].loc[:,[1,8,9,0,6,7,2,3,4,5,10,11]]
    blast_direct_orientation[12] = '+'
    blast_direct_orientation.shape
    blast_reverse_orientation = blast[blast[8]>blast[9]].loc[:,[1,9,8,0,6,7,2,3,4,5,10,11]]
    blast_reverse_orientation[12] = '-'
    blast_reverse_orientation.columns = [1,8,9,0,6,7,2,3,4,5,10,11,12]
    blast_reorder = pd.concat([blast_direct_orientation,blast_reverse_orientation])
    blast_reorder.columns = ['Ref', 'R_Start', 'R_End','Query', 'Q_Start', 'Q_End','%id','alignment_len', 'mismatch', 'gap', 'Evalue', 'score', 'Orientation']
    return blast_reorder

'''Used for remving overlaps'''
def remove_overlap(blast_reorder): 
    RefID = blast_reorder['Ref'].unique().tolist()
    blast_all = pd.DataFrame()
    for ids in RefID:
        #ids= "chr9_modified"
        ids_blast = blast_reorder.loc[blast_reorder['Ref']==ids].sort_values(by='R_Start').reset_index(drop=True) # slice arrays of same ID, sort by start position
        for i in range(1,len(ids_blast.index)): 
            #print(i)
            if int(ids_blast.loc[i-1,'R_End']) < int(ids_blast.loc[i,'R_Start']):
                ids_blast.loc[i-1,'alignment_len_nonoverlap'] = ids_blast.loc[i-1,'alignment_len'] 
                ids_blast.loc[i-1,'R_End_nonoverlap'] = ids_blast.loc[i-1,'R_End'] 
            else:
                ids_blast.loc[i-1,'alignment_len_nonoverlap'] = ids_blast.loc[i-1,'alignment_len'] - abs(int(ids_blast.loc[i-1,'R_End']) - int(ids_blast.loc[i,'R_Start'])) -1
                ids_blast.loc[i-1,'R_End_nonoverlap'] = ids_blast.loc[i-1,'R_End'] - abs(int(ids_blast.loc[i-1,'R_End']) - int(ids_blast.loc[i,'R_Start'])) -1
        ids_blast.loc[len(ids_blast.index)-1,'alignment_len_nonoverlap'] = ids_blast.loc[len(ids_blast.index)-1,'alignment_len']
        ids_blast.loc[len(ids_blast.index)-1,'R_End_nonoverlap'] = ids_blast.loc[len(ids_blast.index)-1,'R_End']
        blast_all = pd.concat([blast_all, ids_blast]).reset_index(drop=True)
    return blast_all

'''Used for clustering and filtering repeats'''    
'''1) cluster repeats by interval spacing 2)filter cluster by length of total repeats in an array and length of an array 3) filter cluster by length of individual repeats in an array'''
'''output repeat pos,size,and structure in each genome, and filtered blast result'''
def cluster_filter(blast_reorder_nonoverlap): 
    RefID=blast_reorder_nonoverlap['Ref'].unique()
    repeat_arrays = pd.DataFrame()
    grouplen_cutoff = 2000 #length of total repeats in an array
    grouplen_total_cutoff = 10000 #length of an array
    space =100000 #interval spacing to define individual repeat array
    individual_repeatlen_cutoff =500 #length of each repeat in an array
    ids_blast_filtered = pd.DataFrame() #blast output filtered according to repeat array length, repeat length in array, and individual repeat length in an array
    for ids in RefID:
        #ids= "Super-Scaffold_1026"
        #cluster
        ids_blast = blast_reorder_nonoverlap.loc[blast_reorder_nonoverlap['Ref']==ids,]
        ids_pos = blast_reorder_nonoverlap.loc[blast_reorder_nonoverlap['Ref']==ids,'R_End_nonoverlap']
        group= []
        breakpoint = 0
        for i in range(len(ids_pos.index)):
            if ids_pos.loc[ids_pos.index[i],] - ids_pos.loc[ids_pos.index[i-1],] > space:       
                group.append(ids_pos.index.tolist()[breakpoint:i])
                breakpoint = i
        group.append(ids_pos.index.tolist()[breakpoint:])
        #filter clustered groups by length of an array and length of total repeats in an array
        group_filtered = []
        group_number = 1
        for i in range(len(group)):
            group_len=ids_blast.loc[min(group[i]):max(group[i])]['alignment_len_nonoverlap'].sum()  
            group_len_total = abs(ids_blast.loc[min(group[i]):max(group[i]),'R_Start'].min()-ids_blast.loc[min(group[i]):max(group[i]),'R_End_nonoverlap'].max())
            if group_len > grouplen_cutoff and group_len_total > grouplen_total_cutoff:
                #print(i)
                
                group_filtered.append(group[i])
                
                repeat_content = [ids,group_number,ids_blast.loc[min(group[i]):max(group[i]),'R_Start'].min(),ids_blast.loc[min(group[i]):max(group[i]),'R_End_nonoverlap'].max(),group_len_total, \
                                                                                               group_len,round(group_len/group_len_total,3)]
                if len(ids_blast.loc[min(group[i]):max(group[i])]['Query'].unique())>1:
                    repeat_content.extend(['Mixed'])
                else:
                    repeat_content.extend(['Single'])
        #further filter clustered groups by individual repeat length in an array
                for x in ids_blast.loc[min(group[i]):max(group[i])]['Query'].unique():
                    for y in ids_blast.loc[min(group[i]):max(group[i])].loc[ids_blast['Query']==x]['Orientation'].unique():
                            if ids_blast.loc[min(group[i]):max(group[i])].loc[(ids_blast['Query']==x) & (ids_blast['Orientation']==y)]['alignment_len_nonoverlap'].sum() > individual_repeatlen_cutoff:
                                repeat_content.extend([x,y,ids_blast.loc[min(group[i]):max(group[i])].loc[(ids_blast['Query']==x) & (ids_blast['Orientation']==y)]['alignment_len_nonoverlap'].sum()])
                                ids_blast_filtered = pd.concat([ids_blast_filtered, ids_blast.loc[min(group[i]):max(group[i])].loc[(ids_blast['Query']==x) & (ids_blast['Orientation']==y)]])
                                
                repeat_arrays = pd.concat([repeat_arrays,pd.Series(repeat_content)], axis=1) 
                group_number += 1
    return repeat_arrays.T,ids_blast_filtered

'''Used for identifying gap distribution and the presence of 100Ns in repeat array'''
def map_gap(gap_file,repeat_sum):
    gap = pd.read_csv(gap_file, delimiter="\t",header=None)
    gap=gap.rename(columns={gap.columns[0]: "Chr", gap.columns[1]: "Start", gap.columns[2]: "End"})
    gap['Size'] = gap['End'] - gap['Start'] 
    #gap.head()
    repeat_sum.insert (5, "100N",pd.Series(["NA"] * repeat_sum.shape[0]))
    repeat_sum.insert (6, "Ngap Size",pd.Series([0] * repeat_sum.shape[0]))
    repeat_sum.insert (7, "Ngap Percentage",pd.Series([0] * repeat_sum.shape[0]))
    for i in repeat_sum.index:
        for j in gap.index:
            if gap.loc[j,"Chr"] == repeat_sum.loc[i,"Chr"] and gap.loc[j,"End"] <= repeat_sum.loc[i,"End"] and gap.loc[j,"Start"] >= repeat_sum.loc[i,"Start"]:
                #print(gap.loc[j,"Size"])
                if gap.loc[j,"Size"] == 100:
                    repeat_sum.loc[i,"100N"] = "Y"
                repeat_sum.loc[i,"Ngap Size"] += gap.loc[j,"Size"]
                repeat_sum.loc[i,"Ngap Percentage"] = round(repeat_sum.loc[i,"Ngap Size"]/repeat_sum.loc[i,"Array size (bp)"],3)
    return repeat_sum


def main(blast_file,gap_file,output_dir):
    blast_reorder = reorder(blast_file)
    blast_reorder_nonoverlap = remove_overlap(blast_reorder)
    repeat_sum,blast_filtered = cluster_filter(blast_reorder_nonoverlap)
    repeat_sum=repeat_sum.rename(columns={repeat_sum.columns[0]: "Chr", repeat_sum.columns[1]: "Group",repeat_sum.columns[2]: "Start", repeat_sum.columns[3]: "End",repeat_sum.columns[4]: "Array size (bp)", repeat_sum.columns[5]: "Repeat content (bp)", \
                           repeat_sum.columns[6]: "Repeat Percentage",repeat_sum.columns[7]: "Repeat Structure", repeat_sum.columns[8]: "Repeat1",repeat_sum.columns[9]: "Orientation1",repeat_sum.columns[10]: "Size1"})
    repeat_sum=repeat_sum.reset_index(drop=True)
    repeat_sum_gap = map_gap(gap_file,repeat_sum)
    os.chdir(output_dir)
    repeat_sum_gap.sort_values(by=['Chr']).to_csv("Repeat_content_sum.csv",na_rep='',index=False)
    for x in blast_reorder_nonoverlap['Query'].unique().tolist():
        x_reverse = blast_reorder_nonoverlap[(blast_reorder_nonoverlap['Query']==x)&(blast_reorder_nonoverlap['Orientation']=='-')].loc[:,['Ref','R_Start','R_End']]
        x_reverse.to_csv("Nonoverlap_{}".format(x) + "_reverse.bed",index=False, sep='\t',header=False)
        x_forward = blast_reorder_nonoverlap[(blast_reorder_nonoverlap['Query']==x)&(blast_reorder_nonoverlap['Orientation']=='+')].loc[:,['Ref','R_Start','R_End']]
        x_forward.to_csv("Nonoverlap_{}".format(x) + "_forward.bed",index=False, sep='\t',header=False)     
    for x in blast_filtered['Query'].unique().tolist():
        x_reverse = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Orientation']=='-')].loc[:,['Ref','R_Start','R_End']]
        x_reverse.to_csv("Filtered_{}".format(x) + "_reverse.bed",index=False, sep='\t',header=False)
        x_forward = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Orientation']=='+')].loc[:,['Ref','R_Start','R_End']]
        x_forward.to_csv("Filtered_{}".format(x) + "_forward.bed",index=False, sep='\t',header=False)     
    repeat_ID = blast_filtered['Query'].unique()
    chr_ID = blast_filtered['Ref'].unique()
    repeat_total_size = []
    for x in repeat_ID:        
        for y in chr_ID: 
            forward_size = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Ref']==y)&(blast_filtered['Orientation']=="+")]["alignment_len_nonoverlap"].sum()
            reverse_size = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Ref']==y)&(blast_filtered['Orientation']=="-")]["alignment_len_nonoverlap"].sum()
            repeat_total_size.append([y,x,"+/-",forward_size+reverse_size,"+",forward_size,"-",reverse_size])
        forward_size_sum = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Orientation']=="+")]["alignment_len_nonoverlap"].sum()
        reverse_size_sum = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Orientation']=="-")]["alignment_len_nonoverlap"].sum()
        repeat_total_size.append(["Genome",x,"+/-",forward_size_sum+reverse_size_sum,"+",forward_size_sum,"-",reverse_size_sum])
        forward_size_chrs = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Ref'].str.contains('chr'))&(blast_filtered['Orientation']=="+")]["alignment_len_nonoverlap"].sum()
        reverse_size_chrs = blast_filtered[(blast_filtered['Query']==x)&(blast_filtered['Ref'].str.contains('chr'))&(blast_filtered['Orientation']=="-")]["alignment_len_nonoverlap"].sum()
        repeat_total_size.append(["Chrs",x,"+/-",forward_size_chrs+reverse_size_chrs,"+",forward_size_chrs,"-",reverse_size_chrs])
    pd.DataFrame(repeat_total_size).sort_values(by=[0,1]).to_csv("Repeat_size_sum.csv",na_rep='',index=False,header=['Chr','Repeat','Total Orientation','Total Size','Orientation1','Size1','Orientation2','Size2'])
    largest_fully_assembled_array = repeat_sum_gap[repeat_sum_gap['100N']!="Y"].sort_values(by=['Array size (bp)'],ascending=False).loc[0]
    largest_fully_assembled_array_withoutNgap = repeat_sum_gap[(repeat_sum_gap['100N']!="Y")&(repeat_sum_gap['Ngap Size']==0)].sort_values(by=['Array size (bp)'],ascending=False).loc[0]
    #print("The largest fully assembled repeat array is " + largest_fully_assembled_array[0]+" located from " + str(int(largest_fully_assembled_array[2])) +"bp to " + str(int(largest_fully_assembled_array[3])) + " bp.")
    #print("The largest fully assembled repeat array without Ngap is " + largest_fully_assembled_array_withoutNgap[0]+ " located from " + str(int(largest_fully_assembled_array_withoutNgap[2])) +"bp to " + str(int(largest_fully_assembled_array_withoutNgap[3])) + " bp.")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Repeat array content analyses')
    parser.add_argument('--blast_file', metavar='', required=True,
                        help='Input blast file (m8 format)')
    parser.add_argument('--gap_file', metavar='', required=True,
                        help='Gaps in genome (bed file form sat)')
    parser.add_argument('--output_dir', metavar='', required=True,
                        help='Output directory')
    args = parser.parse_args()
    main(blast_file=args.blast_file, gap_file=args.gap_file, output_dir=args.output_dir)
    
