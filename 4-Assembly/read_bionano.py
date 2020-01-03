import pandas as pd 
import numpy as np

'''Read BioNano files'''
def read_smap(smap_file):
    smap = pd.read_csv(smap_file, delimiter="\t",comment="#",header=None)
    smap = smap.loc[:,[1,4,5,2,3,6,7,10,11,13,14,15,16,9,24,26,17,8]]
    smap.columns = ['Query_ID', 'Q_Start', 'Q_End','Ref_ID1', 'Ref_ID2','R_Start', 'R_End','XmapID1','XmapID2', 'QryStartIdx', 'QryEndIdx', 'RefStartIdx', 'RefEndIdx','SV_Type', 'SV_Size','Orientation', 'Zygosity', 'Confidence']
    smap_pav = smap[(smap['Zygosity']=="homozygous") & (smap['Confidence']>=0.1) & (smap['SV_Size']>=1000) & (smap['SV_Type'].str.contains("insertion|deletion|duplication")==True)]
    smap_trans = smap[(smap['Zygosity']=="homozygous") & (smap['SV_Type'].str.contains("trans")==True)]
    return smap_pav.reset_index(drop=True), smap_trans.reset_index(drop=True),smap

def read_xmap(xmap_file):
    xmap = pd.read_csv(xmap_file, delimiter="\t",comment="#",header=None).loc[:,[0,1,3,4,2,5,6,7,10,11,13]]
    xmap.columns = ['CMapId','Query_ID', 'Q_Start', 'Q_End','Ref_ID1', 'R_Start', 'R_End','Orientation','Q_Len', 'R_Len', 'Alignment_String']
    return xmap.reset_index(drop=True)

def read_cmap(cmap_file):
    cmap = pd.read_csv(cmap_file,delimiter="\t",comment="#",header=None)
    cmap.columns = ['CMapId','ContigLength', 'NumSites', 'SiteID','LabelChannel', 'Position','StdDev','Coverage','Occurrence']
    return cmap.reset_index(drop=True)

def read_key(key_file):
    key = pd.read_csv(key_file,delimiter="\t",comment="#",header=0)
    return key.reset_index(drop=True)

def read_cut(cut_file,keys):
    cut = pd.read_csv(cut_file, delimiter="\t",comment="#",header=0)
    cut.insert(6, "Sequence_ID", 'NA') 
    # map cut file ID to sequence ID, with Key file
    for i in cut.index:
        key_id = cut.loc[i,'oldId']
        cut.loc[i,'Sequence_ID'] = keys[keys['CompntId']==key_id]['CompntName'].values[0]
    return cut
