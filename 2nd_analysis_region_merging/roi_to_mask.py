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
###########





index_csf = legend [ 'Structure' ] == 'Amygdala' 
index_csf+= legend [ 'Structure' ] == 'Basal Lateral Amygdala' 
index_csf+= legend [ 'Structure' ] == 'Posterolateral_Cortical_Amygdaloid_Area' 



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
nib.save(file_result,mypath + '2nd_analysis_region_merging/Amygdala.nii.gz'  )

#label_path_res= mypath+'chass_symmetric3_labels_PLI_res.nii.gz'
#os.system('/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i ' +mypath + 'olfac_amyg_hippo.nii.gz -n NearestNeighbor -r '+label_path_res+" -o "+mypath + 'olfac_amyg_hippo_0P3.nii.gz') 




####### remaining 

index_csf = legend [ 'Structure' ] != '' 


values= [
         'Amygdala',
         'Posterolateral_Cortical_Amygdaloid_Area', 
         'Basal Lateral Amygdala',
         'Cerebellar_Cortex',
         'Dentate_(Lateral)_Nucleus_of_Cerebellum',
         'Interposed_Nucleus_of_Cerebellum',
         'Fastigial_Medial_Dorsolateral_Nucleus_of_Cerebellum',
         'Fastigial_Medial_Nucleus_of_Cerebellum',
         'Thalamus_Rest',
'Ventral_Thalamic_Nuclei',
'Latero_Dorsal_Nucleus_of_Thalamus',
'Medial_Geniculate_Nucleus',
'Lateral_Geniculate_Nucleus',
'Zona_Incerta',
'Reticular_Nucleus_of_Thalamus',
'Ventral_Claustrum',
        'Claustrum',
'Dorsal_Claustrum', 
'Caudomedial_Entorhinal_Cortex',
'Dorsal_Intermediate_Entorhinal_Cortex',
'Dorsolateral_Entorhinal_Cortex',
'Medial_Entorhinal_Cortex',
'Ventral_Intermediate_Entorhinal_Cortex',
'Secondary_Visual_CortexLateral_Area',
'Secondary_Visual_Cortex_Mediolateral_Area',
'Secondary_Visual_Cortex_Mediomedial_Area',
'Primary_Visual_Cortex',
'Primary_Visual_Cortex_Binocular_Area',
'Primary_Visual_Cortex_Monocular_Area',
'Primary_Somatosensory_Cortex',
'Primary_Somatosensory_Cortex_Barrel_Field',
'Primary_Somatosensory_Cortex_Dysgranular_Zone',
'Primary_Somatosensory_CortexForelimb_Region',
'Primary_Somatosensory_Cortex_Hindlimb_Region',
'Primary_Somatosensory_Cortex_Jaw_Region',
'Primary_Somatosensory_Cortex_Shoulder_Region',
'Primary_Somatosensory_Cortex_Trunk_Region',
'Primary_Somatosensory_Cortex_Upper_Lip_Region',
'Primary_Motor_Cortex',
'Secondary_Motor_Cortex',
'Secondary_Auditory_Cortex_Dorsal_Part',
'Secondary_Auditory_Cortex_Ventral_Part',
'Cingulate_Cortex_Area_24a',
'Cingulate_Cortex_Area_24a_prime',
'Cingulate_Cortex_Area_24b',
'Cingulate_Cortex_Area_24b_prime',
'Cingulate_Cortex_Area_29a',
'Cingulate_Cortex_Area_29b',
'Cingulate_Cortex_Area_29c',
'Cingulate_Cortex_Area_30',
'Cingulate_Cortex_Area_32'
         ]



index_csf = np.transpose([~legend['Structure'].isin(values)])

vol_index_csf = legend[index_csf]
vol_index_csf = vol_index_csf['index2']







label_nii_csf_data =label_nii.get_fdata()*0

for csf in vol_index_csf:
    #print(csf)
    label_nii_csf_data[  data_label == int(csf)] = 1
    
    
    
file_result= nib.Nifti1Image(label_nii_csf_data, label_nii.affine, label_nii.header)
nib.save(file_result,mypath + '2nd_analysis_region_merging/remaining.nii.gz'  )







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