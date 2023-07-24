#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:37:17 2023

@author: ali
"""

import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd
import openpyxl

mypath= '/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/'

#construct csf and wm mask
label_path= mypath + 'chass_symmetric3_labels_PLI_res.nii.gz'
label_nii=nib.load(label_path)
label_nii.shape
data_label=label_nii.get_fdata()
roi_list=np.unique(data_label)
roi_list = roi_list[1:]

proportion = []

for roi in roi_list:
    np.mean(roi==data_label)
    proportion.append(    np.mean(roi==data_label) )

top_number_regions= 20 
top_indecies= sorted(range(len(proportion)), key=lambda i: proportion[i], reverse=True)[:20]
# these can be only usefull in 324 Atlas file


path_atlas_legend ="/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/new_atlas.xlsx"
legend  = pd.read_excel(path_atlas_legend)
new_atlas_index_top = [int(i)-1 for i in top_indecies]
atlas_top_indecies = legend.loc[new_atlas_index_top ]
atlas_top_indecies.to_excel("/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/top_regions.xlsx")

''' here the all regions individual files
path_atlas_legend ="/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/CHASSSYMM3AtlasLegends.xlsx"
legend  = pd.read_excel(path_atlas_legend)

set(legend['index2'])  & set(roi_list)
index = [i for i, item in enumerate(legend['index2']) if item in set(roi_list)]
new_leg = legend.iloc[index]

new_leg.to_excel( "/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/CHASSSYMM3AtlasLegends_324.xlsx" )


List_of_ROIs= np.unique(legend [ 'Structure' ] )
ind_total = List_of_ROIs == 'TotalBrain'
List_of_ROIs =  np.delete(List_of_ROIs , ind_total)

for reg in List_of_ROIs:

    index_csf = legend [ 'Structure' ] == reg 
#index_csf+= legend [ 'Structure' ] == 'Hippocampus' 
#index_csf+= legend [ 'Structure' ] == 'Hippocampus' 
#index_csf+= legend [ 'Structure' ] == 'Lateral_Olfactory_Tract' 


#index_wm = legend [ 'Structure' ] == 'Hippocampus'

    vol_index_csf = legend[index_csf]
    vol_index_csf = vol_index_csf['index2']



#vol_index_wm  = legend[index_wm]
#vol_index_wm  = vol_index_wm['index2']




    label_nii_csf_data =label_nii.get_fdata()*0

    for csf in vol_index_csf:
    #print(csf)
        label_nii_csf_data[  data_label == int(csf)] = 1
    
    
    
    file_result= nib.Nifti1Image(label_nii_csf_data, label_nii.affine, label_nii.header)
    nib.save(file_result,mypath +'each_region/'+ reg+'.nii.gz'  )

#label_path_res= mypath+'chass_symmetric3_labels_PLI_res.nii.gz'
#os.system('/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i ' +mypath + 'olfac_amyg_hippo.nii.gz -n NearestNeighbor -r '+label_path_res+" -o "+mypath + 'olfac_amyg_hippo_0P3.nii.gz') 




'''

label_nii_wm_data =label_nii.get_fdata()*0

for wm in vol_index_wm:
    #print(csf)
    label_nii_wm_data[  data_label == int(wm)] = 1
    
    
    
file_result= nib.Nifti1Image(label_nii_wm_data, label_nii.affine, label_nii.header)
nib.save(file_result,mypath + 'chass_symmetric3/wm_mask.nii.gz'  )
os.system('/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i ' +mypath +'chass_symmetric3/wm_mask.nii.gz -n NearestNeighbor -r '+label_path_res+" -o "+mypath +'chass_symmetric3/wm_mask_0p3.nii.gz') 


########### making a mask out of labels

label_path= mypath +'/chass_symmetric3/chass_symmetric3_labels_PLI_res.nii.gz'
label_nii=nib.load(label_path)
mask_labels_data = label_nii.get_fdata()
mask_labels = np.unique(mask_labels_data)
mask_labels=np.delete(mask_labels, 0)
mask_of_label =label_nii.get_fdata()*0

for vol in mask_labels:
    mask_of_label[  mask_labels_data == int(vol)] = 1
    
file_result= nib.Nifti1Image(mask_of_label, label_nii.affine, label_nii.header)
nib.save(file_result,mypath + 'chass_symmetric3/mask_of_label.nii.gz'  )    
'''   