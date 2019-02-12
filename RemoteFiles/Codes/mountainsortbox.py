import glob, os, shutil, re
import numpy as np

def prepare_mountainsort_files(mdafiles=(), rawfolder='/tmp/.x2go-jsee/media/disk/_cygdrive_C_Users_Jermyn_BANSHE1/MSInput',
                               datapath='/home/jsee/MountainSort/datasets/input', outputpath='/home/jsee/MountainSort/datasets/output',
                               prvopt=False, combineopt=False, combinelist='auto'):

    print('Creating input folders...')
    fplinput, fploutput, mdafiles = create_dataset_folders(mdafiles, rawfolder, datapath, outputpath)

    c = 1
    for f,m in zip(fplinput, mdafiles):
        print('Copying files for dataset {} of {}...'.format(c, len(mdafiles)))
        copy_input_files(f, m, prvopt)
        c += 1

    if prvopt:
        print('Writing bash script...')
        create_bash_script(fplinput, fploutput, mdafiles)

    if combineopt:

        if combinelist == 'auto':
            groupedinput, groupedoutput, siteindex = auto_group_same_sites(fplinput, fploutput)
        else:
            pass    # will code this when need arises

        print('Creating combined folders...')
        combfolders = create_combined_output_folder(groupedoutput)

        c = 1
        for si, cf in zip(siteindex, combfolders):
            print('Writing combined raw .mda files for dataset {} of {}...'.format(c, len(combfolders)))
            mdalist = [mdafiles[i] for i in si]
            combine_raw_mda_files(mdalist, cf)
            c += 1

        return groupedinput, groupedoutput, combfolders

    else:

        return fplinput, fploutput



def create_dataset_folders(mdafiles=(), rawfolder='/tmp/.x2go-jsee/media/disk/_cygdrive_C_Users_Jermyn_BANSHE1/MSInput',
                               datapath='/home/jsee/MountainSort/datasets/input', outputpath='/home/jsee/MountainSort/datasets/output'):

    # get all .mda files in shared folder if no mdafile input
    if not any(mdafiles):
        mdafiles = glob.glob(os.path.join(rawfolder, '*.mda'))
    else:
        mdafiles = [os.path.join(rawfolder, m) for m in mdafiles]

    # datapath = '/home/jsee/MountainSort/datasets/input'
    # outputpath = '/home/jsee/MountainSort/datasets/output'
    fplinput = []
    fploutput = []

    for i in mdafiles:

        # check for subfolder in dataset; if it doesn't exist, make the folder
        basefilename = os.path.basename(i).split('.')[0]
        fullinputpath = os.path.join(datapath, basefilename)

        outputfolder = os.path.basename(re.sub('mountainsort', 'firings', i)).split('.')[0]
        fulloutputpath = os.path.join(outputpath, outputfolder)

        if not os.path.isdir(fullinputpath):
            os.mkdir(fullinputpath)

        if not os.path.isdir(fulloutputpath):
            os.mkdir(fulloutputpath)

        fplinput.append(fullinputpath)
        fploutput.append(fulloutputpath)

    return fplinput, fploutput, mdafiles


def copy_input_files(folderpath, mdafile, prvopt):


    prbfolder = '/home/jsee/MountainSort/prbfiles'

    prbname = re.search('(?<=min-)\S+(?=-fs)', folderpath).group(0)
    prbfile = glob.glob(os.path.join(prbfolder, prbname + '.csv'))
    assert len(prbfile) == 1

    shutil.copyfile(prbfile[0], os.path.join(folderpath, 'geom.csv'))

    # paramsfile = os.path.join(prbfolder, 'params.json')
    # shutil.copyfile(paramsfile, os.path.join(folderpath, 'params.json'))

    if not prvopt:
        if not os.path.isfile(os.path.join(folderpath, 'raw.mda')):
            shutil.copyfile(mdafile, os.path.join(folderpath, 'raw.mda'))
        else:
            print('raw.mda file in {} already copied! Skipping...'.format(folderpath))


def combine_raw_mda_files(filelist, outfolder):

    if os.path.isfile(outfolder + '/combined_raw.mda'):
        print('{} already exists! Skipping...'.format(outfolder+'/combined_raw.mda'))
    else:
        paramlist = []
        nsamples = 100000000

        for file in filelist:
            with open(file, 'r') as fid:
                paramlist.append(np.fromfile(fid, 'int32', 5))

        paramzip = list(zip(*paramlist))

        # Ensure all parameters are the same (except for length, which is the last param)
        assert(all([len(list(set(i))) == 1 for i in paramzip[0:-1]]))

        newlength = np.sum(paramzip[-1])
        newparam = [p[0] for p in paramzip]
        newparam[-1] = newlength
        newparam = np.array(newparam, dtype='int32')

        with open(outfolder + '/combined_raw.mda', 'w') as fidout:
            newparam.tofile(fidout)

            for file in filelist:

                with open(file, 'r') as fidin:
                    __ = np.fromfile(fidin, 'int32', 5)

                    c = 1

                    while True:

                        towrite = np.fromfile(fidin, 'int16', nsamples)

                        if towrite.size == 0:
                            break

                        towrite.tofile(fidout)
                        print('{} samples processed'.format(c * nsamples))
                        c += 1


def auto_group_same_sites(fplinput, fploutput):

    sitenum = [int(re.search('(?<=site)\d{1,2}', i).group(0)) for i in fploutput]
    uniquesites = list(set(sitenum))
    uniquesites.sort()

    groupedinput = []
    groupedoutput = []
    siteindex = []
    for s in uniquesites:

        groupedinput.append([i for i, site in zip(fplinput, sitenum) if site == s])
        groupedoutput.append([i for i, site in zip(fploutput, sitenum) if site == s])
        siteindex.append([i for i, site in enumerate(sitenum) if site == s])

    return groupedinput, groupedoutput, siteindex

def create_combined_output_folder(groupedsites):

    combfolders = []

    for i in groupedsites:

        combstim = [re.search('(?<=db-)\w+(?=-\d{1,3}min)', j).group(0) for j in i]
        combstim = '_'.join(combstim)
        stimlen = str(sum([int(re.search('(?<=-)\d+(?=min)', j).group(0)) for j in i])) + 'min'
        replace = combstim + '-' + stimlen

        sampfolder = i[0]
        toreplace = re.search('(?<=db-)\S+min(?=-)', sampfolder).group(0)
        groupedfolder = re.sub(toreplace, replace, sampfolder)

        combfolders.append(groupedfolder)

        if not os.path.isdir(groupedfolder):
            os.mkdir(groupedfolder)

    return combfolders


def copy_output_to_shared_folder(sharedfolder, outputfolders, filenames=()):

    if not filenames:
        filenames = ('cluster_metrics.json', 'firings.mda', 'firings_uncurated.mda')

    print('Transferring output files...')

    for of in outputfolders:

        basefolder = os.path.basename(of)
        fullopfolder = os.path.join(sharedfolder, basefolder)

        if not os.path.isdir(fullopfolder):
            os.mkdir(fullopfolder)

        [shutil.copyfile(os.path.join(of, fn), os.path.join(fullopfolder, fn)) for fn in filenames]

# def add_datafolder_to_datasetstxt(folderpaths):
#
#     # add '/' to each folderpath if it does not end with a '/'
#     # folderpaths = [i + '/' if i[-1] != '/' else i for i in folderpaths]
#
#     os.chdir('/home/jsee/mountainlab/data')
#
#     sortnames = []
#
#     if os.path.isfile('datasets.txt'):
#
#         with open('datasets.txt', 'r') as f:
#             dstext = f.readlines()
#
#         for i in folderpaths:
#
#             relpath = re.split('(?<=data)/(?=datasets)', i)[-1]
#             stim = re.search('(?<=db-)\w+(?=-)', relpath).group(0)
#             base = re.search('(?<=datasets/)\S+site\d{1,2}(?=-)', relpath).group(0)
#             name = base + '-' + stim
#             sortnames.append(name)
#
#             if not any([relpath in j for j in dstext]):
#
#                 addtext = name + ' ' + relpath
#                 dstext.append(addtext + '\n')
#
#         # delete old textfile
#         os.remove('datasets.txt')
#
#     else:
#
#         dstext = []
#
#         for i in folderpaths:
#             relpath = re.split('(?<=data)/(?=datasets)', i)[-1]
#             stim = re.search('(?<=db-)\w+(?=-)', relpath).group(0)
#             base = re.search('(?<=datasets/)\S+site\d{1,2}(?=-)', relpath).group(0)
#             name = base + '-' + stim
#             addtext = name + ' ' + relpath
#             dstext.append(addtext + '\n')
#             sortnames.append(name)
#
#
#     # write new textfile
#     with open('datasets.txt', 'w') as f:
#         f.writelines(dstext)
#
#     return sortnames


def create_bash_script(fplinput, fploutput, mdafiles, freqmin=300, freqmax=6000, adjradius=200):

    assert len(fplinput) == len(mdafiles)
    assert len(fploutput) == len(mdafiles)

    with open('/home/jsee/MountainSort/datasets/input/batchMS.sh', 'w') as f:
        f.write('#!/bin/bash\n\n')
        # f.write('cd /home/jsee/MountainSort/datasets/input\n')
        f.write('set -e\n\n')


        for i in range(len(fplinput)):

            # rawfile = os.path.join(rawpath, mdafiles[i])
            fs = int(re.search('(?<=fs)\d{4,5}', mdafiles[i]).group(0))
            prvfile = os.path.join(fplinput[i], 'raw.mda.prv')
            geomfile = os.path.join(fplinput[i], 'geom.csv')
            filtfile = os.path.join(fploutput[i], 'filt.mda.prv')
            whitenfile = os.path.join(fploutput[i], 'pre.mda.prv')
            templatefile = os.path.join(fploutput[i], 'templates.mda')
            outputfile = os.path.join(fploutput[i], 'firings.mda')
            metricsfile1 = os.path.join(fploutput[i], 'metrics1.json')
            metricsfile2 = os.path.join(fploutput[i], 'metrics2.json')
            metricsfile3 = os.path.join(fploutput[i], 'metrics3.json')
            pairmetricsfile = os.path.join(fploutput[i], 'pair_metrics.json')
            clusterfile = os.path.join(fploutput[i], 'cluster_metrics.json')

            # create prv file
            f.write('ml-prv-create {} {}\n\n'.format(mdafiles[i], prvfile))

            # bandpass filter
            f.write('ml-run-process ephys.bandpass_filter \\\n\t--inputs timeseries:{} \\\n\t--outputs timeseries_out:{}'
                    ' \\\n\t--parameters samplerate:{} freq_min:{} freq_max:{}\n\n'
                    .format(prvfile, filtfile, fs, freqmin, freqmax))

            # whiten data
            f.write('ml-run-process ephys.whiten \\\n\t--inputs timeseries:{} \\\n\t--outputs timeseries_out:{}\n\n'
                    .format(filtfile, whitenfile))

            # run sort
            f.write('ml-run-process ms4alg.sort \\\n\t--inputs timeseries:{} geom:{} \\\n\t--outputs firings_out:{} '
                    '\\\n\t--parameters detect_sign:-1 adjacency_radius:{} detect_threshold:3\n\n'
                    .format(whitenfile, geomfile, outputfile, adjradius))

            # compute metrics
            f.write('ml-run-process ms3.cluster_metrics \\\n\t--inputs timeseries:{} firings:{} \\\n\t--outputs '
                    'cluster_metrics_out:{} \\\n\t--parameters samplerate:{}\n\n'
                    .format(prvfile, outputfile, metricsfile1, fs))

            f.write('ml-run-process ms3.isolation_metrics \\\n\t--inputs timeseries:{} firings:{} \\\n\t--outputs '
                    'metrics_out:{} pair_metrics_out:{} \\\n\t--parameters compute_bursting_parents:true\n\n'
                    .format(prvfile, outputfile, metricsfile2, pairmetricsfile))

            # f.write('ml-run-process ephys.compute_cluster_metrics \\\n\t--inputs firings:{} timeseries:{} \\\n\t'
            #         '--outputs metrics_out:{} \\\n\t--parameters samplerate:{} refrac_msec:1.5\n\n'
            #         .format(outputfile, prvfile, metricsfile3, fs))

            f.write('ml-run-process ms3.combine_cluster_metrics \\\n\t--inputs metrics_list:[{}, {}, {}] \\\n\t'
                    '--outputs metrics_out:{}\n\n'
                    .format(metricsfile1, metricsfile2, metricsfile3, clusterfile))

            # get waveforms
            f.write('ml-run-process ephys.compute_templates \\\n\t--inputs timeseries:{} firings:{} \\\n\t--outputs '
                    'templates_out:{} \\\n\t--parameters clip_size:150\n\n\n\n'
                    .format(prvfile, outputfile, templatefile))



            f.write('ml-run mountainsort3.mlp sort --raw={} --geom={} --firings_out={} --firings_original_out={} --cluster_metrics_out={} --_params=params.json --curate=true\n'.format(prvfile, geomfile, outputfile, outputfileori, clusterfile))
