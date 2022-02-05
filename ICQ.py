#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import numpy as np
import imageio
import pandas as pd
import os
import fnmatch
import seaborn as sns
import datetime

#%%
fpath = 'G:/My Drive/ECM manuscript/github codes/ICQ/sample_data/input_files/'
dfpath = 'G:/My Drive/ECM manuscript/github codes/ICQ/sample_data/output_files/'
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')
toa = str(datetime.datetime.today()).split()
today = toa[0]
now = toa[1]
timestamp = today.replace('-','')+'-'+now.replace(':','')[:6]

strain_key=pd.DataFrame({('x42', 'GN902', 'MEC-4 LAM-1', 'WT', 0),
                         ('x68', 'GN941', 'MEC-4 NID-1', 'WT',0),
                         ('x69', 'GN940', 'LAM-2 NID-1', 'WT',0),
                         ('x70', 'GN939', 'LAM-2 LAM-1', 'WT',0),
                         ('x81', 'GN966', 'MEC-4 PAT-2', 'WT',0),
                         ('x112', '', 'NID-1 NID-1', 'WT',0),
                         ('x84', 'GN1000', 'MEC-4 LAM-1', 'mec-1(e1526)', 0),
                         ('x83', 'GN999', 'MEC-4 LAM-1', 'nid-1(cg119)', 0),
                         ('x119', 'GN1004', 'MEC-4 NID-1', 'mec-1(e1526)', 0),
                         ('x120', 'GN1005', 'LAM-2 NID-1', 'mec-1(e1526)', 0),
                         ('inj61', 'GN1062', 'MEC-4 NID-1', 'mec-1(pg154)', 0),
                         ('inj64', 'GN1073', 'MEC-4 NID-1', 'mec-1(pg164)', 0)
                         }, columns=['Strain_code', 'Strain','G/R', 'Background', 'n'])
# strain_key=strain_key.set_index('G/R')

mu_per_px = 0.126     #pixels to microns conversion factor

def corr_coeffs(rawG, rawR):
    nG = rawG[7:13]
    nR = rawR[7:13]
    bgG= np.concatenate((rawG[0:6, 0:], rawG[14: , 0:]))
    bgR= np.concatenate((rawR[0:6, 0:], rawR[14: , 0:]))
    G = np.subtract(nG, np.mean(bgG,axis=0))
    R = np.subtract(nR, np.mean(bgR,axis=0))
    dfmG = np.subtract(G,np.mean(G))
    dfmR = np.subtract(R,np.mean(R))
    pdm=np.multiply(dfmG,dfmR)                  #product matrix
    totalpx = np.shape(pdm)[0]*np.shape(pdm)[1] #total number of pixels in product matrix
    pospx = np.count_nonzero(pdm>0)             #total number of positive pixels
    icq=float(pospx)/totalpx - 0.5
    pearson = np.sum(pdm)/np.sqrt(np.sum(np.power(dfmG,2))*np.sum(np.power(dfmR,2)))
    return(icq, pearson)

proximal=[0.05, 0.30]
distal=[0.65, 0.90]


cols_ICQ = ['Strain', 'G/R', 'Background', 'ImageID', 'Length', 'Segment', 'ICQ', 'Pearson coefficient']
df_ICQ = pd.DataFrame()

#%%
for img in imgfiles:
    
    ID = img.split('_')[1].split('-')[0]
    row_index=strain_key[(strain_key['Strain']==ID)|(strain_key['Strain_code']==ID)].index[0]
    strain = strain_key.loc[row_index,'Strain']
    proteinpair = strain_key.loc[row_index,'G/R']
    bg = strain_key.loc[row_index,'Background']
    n = strain_key.loc[row_index,'n'] + 1
    strain_key.at[row_index,'n']=n

    imgG = imageio.imread(fpath+img)[:,:,1]           #import image and store it in a list of lists
    imgR = imageio.imread(fpath+img)[:,:,0]           #import image and store it in a list of lists
        
    imsize = np.shape(imgG)                  #calculate image size
    dist = np.arange(0,imsize[1])
    d=dist*mu_per_px
    norm_dist = np.divide(dist,imsize[1]-1)

    icq_full, pearson_full = corr_coeffs(imgG, imgR)
    
    frame1 = pd.DataFrame([[strain, proteinpair, bg, img, imsize[1], 'full', icq_full, pearson_full]], columns=cols_ICQ)
    df_ICQ = df_ICQ.append(frame1)

#%%
#save dataframes
wb = pd.ExcelWriter(dfpath+timestamp+'_Analysis.xlsx', engine='xlsxwriter')
df_ICQ.to_excel(wb, sheet_name='ICQ')
strain_key.to_excel(wb, sheet_name='Summarized strain info')
wb.save()
df_ICQ.to_pickle(dfpath+timestamp+'_ICQ.pkl')
   


