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



path_atlas_legend ="/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/CHASSSYMM3AtlasLegends.xlsx"
legend  = pd.read_excel(path_atlas_legend)

set(legend['index2'])  & set(roi_list)
index = [i for i, item in enumerate(legend['index2']) if item in set(roi_list)]
new_leg = legend.iloc[index]

new_leg.to_excel( "/Users/ali/Desktop/Jul23/atlas_maker_jayvick_mouse_network/CHASSSYMM3AtlasLegends_324.xlsx" )




index_csf = legend [ 'Structure' ] == 'Amygdala' 
index_csf+= legend [ 'Structure' ] == 'Hippocampus' 
index_csf+= legend [ 'Structure' ] == 'Hippocampus' 
index_csf+= legend [ 'Structure' ] == 'Lateral_Olfactory_Tract' 


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
nib.save(file_result,mypath + 'olfac_amyg_hippo.nii.gz'  )

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