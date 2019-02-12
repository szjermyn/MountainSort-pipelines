import spikesortbox as ssb

outfolder = r'I:\SM_data\nhp_su_anesthetized\ContraIpsiThresh'
folder = r'C:\Users\Jermyn\BansheeShared\MSResults\ContraIpsiThresh'
triggerfolder = r'L:\Anesthetized SM Data'
ssb.batch_save_mountainsort_results(outfolder, folder=folder, triggerfolder=triggerfolder, combineopt=False)