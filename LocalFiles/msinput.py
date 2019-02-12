import spikesortbox as ssb
import os, glob, re


# This section is for preparing mda files for sorting
# Change these:
folderpath = r'L:\Anesthetized SM Data' #folder with subfolders for where .raw files are
subfolderid = '*db_dmr_contra_ipsi_H*' #process subfolders containing this substring. usually a stimulus string.
outfolder = r'C:\Users\Jermyn\BansheeShared\MSInput\Demo'

stim = ()

os.chdir(folderpath)
folders = glob.glob(subfolderid)
folders = [os.path.join(folderpath, i) for i in folders]

ssb.batch_raw2mda(folders, stim, outfolder)



# '''For combinations'''
# folderpath = r'L:\Anesthetized SM Data' #folder with subfolders for where .raw files are
# subfolderid1 = '*dmr_contra_H*' #process subfolders containing this only. usually a stimulus string.
# subfolderid2 = '*dmr_contra_ipsi_H*'
# outfolder = r'C:\Users\Jermyn\BansheeShared\MSInput\combdebug'
#
# stim = ()
#
# os.chdir(folderpath)
# folders1 = glob.glob(subfolderid1)
# del folders1[4]
# sites1 = [re.search('site\d{1,2}_',i).group(0) for i in folders1]
#
# folders2 = glob.glob(subfolderid2)
# sites2 = [re.search('site\d{1,2}_',i).group(0) for i in folders2]
#
# folders2 = [folders2[i] for i in range(len(sites2)) if sites2[i] in sites1]
#
# folders = [[os.path.join(folderpath, i), os.path.join(folderpath, j)] for i,j in zip(folders1, folders2)]
# # folders = folders[:3]
#
# folders = [[r'L:\Anesthetized SM Data\site2_1380um_hp20db_art32db_dmr_contra_H2_2x32_180802_234840',r'L:\Anesthetized SM Data\site2_1380um_hp20db_art32db_dmr_contra_ipsi_H2_2x32_180803_001032']]
# ssb.batch_raw2mda(folders, stim, outfolder, combineopt=True)