import os, glob, re, shutil, json, csv, warnings, itertools
import numpy as np
import scipy.io as sio

def batch_raw2mda(folder, stim=(), outfilefolder=(), splitopt=False, combineopt=False):

    if combineopt and splitopt:
        raise Exception('Only one of combineopt and splitopt can be True')

    for i in folder:
        # os.chdir(i)
        print('Processing {}...'.format(i))
        raw2mda(i, stim, outfilefolder, splitopt, combineopt)


def raw2mda(folder, stim=(), outfilefolder=(), splitopt=False, combineopt=False):

    # combineopt assumes that there is only one stim per folder

    if combineopt and splitopt:
        raise Exception('Only one of combineopt and splitopt can be True')

    if combineopt:
        assert(type(folder) is list)
    else:
        assert(type(folder) is str)

    if not stim:
        if combineopt:
            filesearch = list(itertools.chain.from_iterable([glob.glob(os.path.join(i, '*min*.raw')) for i in folder]))
            temp = [re.search('(?<=db-)\w+(?=-)', i) for i in filesearch]
            stim = [list(dict.fromkeys([i.group(0) for i in temp if i is not None]))]
        else:
            filesearch = glob.glob(os.path.join(folder, '*min*.raw'))
            temp = [re.search('(?<=db-)\w+(?=-)', i) for i in filesearch]
            stim = list(dict.fromkeys([i.group(0) for i in temp if i is not None]))

    elif isinstance(stim, str):
        stim = [stim]

    num_dims = 2
    nsamples = 1000000

    if not outfilefolder:
        outfilefolder = r'C:\Users\Jermyn\BansheeShared\MSInput'

    for i in stim:

        # get file names
        if combineopt:
            files = list(itertools.chain.from_iterable([glob.glob(os.path.join(f, '*-{}-*min*.raw'.format(j))) for j in i for f in folder]))
        else:
            files = glob.glob(os.path.join(folder, '*-{}-*min*.raw'.format(i)))

        # remove trig file
        trigfiles = [f for f in files if 'ADC' in f]
        [files.remove(tf) for tf in trigfiles]

        # check for multiple probes
        usedamps = [re.search('\w(?=-\d{3}.raw)', i).group(0) for i in files]
        uniqamps = set(usedamps)

        #put probes in separate .mda files if required
        if len(uniqamps) > 1:

            if splitopt:
                originaldepth = re.search('\w+(?=um)', files[0]).group(0)
                originalprobe = re.search('(?<=min-)\S+(?=-fs)', files[0]).group(0)

            batchedfiles = []
            outfile = []

            for j in uniqamps:
                index = [jj for jj, kk in enumerate(usedamps) if kk == j]
                batchedfiles.append([files[k] for k in index])
                ampidx = list(uniqamps).index(j)

                if splitopt:

                    depth = re.findall('\w\d{2,5}', originaldepth)[ampidx]
                    probe = re.findall('\w+-\dx\d{1,3}', originalprobe)[ampidx]
                    filebase = re.match('^\S+(?=-[ABCD]-\d{3}.raw)', files[0]).group(0)
                    filebase = re.sub(originaldepth, depth, filebase)
                    filebase = re.sub(originalprobe, probe, filebase)

                else:

                    filebase = re.match('^\S+(?=-[ABCD]-\d{3}.raw)', os.path.basename(files[index[0]])).group(0)

                outfile.append(r'{}\{}-mountainsort.mda'.format(outfilefolder, filebase))


        elif combineopt:

            filestim = [re.search('(?<=db-)\w+(?=-\d{1,2}min)', j).group(0) for j in files]
            batchedfiles = [[files[idx] for idx,j in enumerate(filestim) if j == s] for s in i]
            filebase = re.match('(^\S+db-)\w+(-\d{1,2}min\S+(?=-[ABCD]-\d{3}.raw))', os.path.basename(batchedfiles[0][0]))
            mdastim = '_'.join(i)
            outfile = '{}\\{}{}{}-mountainsort.mda'.format(outfilefolder, filebase.group(1), mdastim, filebase.group(2))


        else:
            # get outfile name
            batchedfiles = [files]
            filebase = re.match('^\S+(?=-[ABCD]-\d{3}.raw)', os.path.basename(batchedfiles[0][0]))
            outfile = ['{}\\{}-mountainsort.mda'.format(outfilefolder, filebase.group(0))]


        if combineopt:

            filelength = sum([int(os.path.getsize(j[0]) / 2) for j in batchedfiles])

            with open(outfile, 'w') as fidout:
                iniparams = np.empty([5, 1], dtype='int32')
                iniparams[0] = -4  # specify int16 type for data
                iniparams[1] = 2  # specify number of bytes for each int16 variable
                iniparams[2] = num_dims  # specify number of dimensions for data
                iniparams[3] = len(batchedfiles[0])
                iniparams[4] = filelength
                iniparams.tofile(fidout)

                for j in batchedfiles:

                    fidin = [open(k, 'r') for k in j]

                    c = 1
                    br = 0
                    while True:
                        towrite = []
                        for k in fidin:
                            towrite.append(np.fromfile(k, 'int16', nsamples))
                            if towrite[0].size == 0:
                                br = 1
                                break
                        if br:
                            break
                        towrite = np.asarray(towrite, 'int16')
                        towrite = towrite.flatten('F')
                        towrite.tofile(fidout)
                        print('{} samples processed...'.format(c * nsamples))
                        c += 1


        else:

            for j, of in zip(batchedfiles,outfile):

                filelength = int(os.path.getsize(j[0])/2)

                with open(of, 'w') as fidout:
                    iniparams = np.empty([5,1], dtype='int32')
                    iniparams[0] = -4 #specify int16 type for data
                    iniparams[1] = 2 #specify number of bytes for each int16 variable
                    iniparams[2] = num_dims #specify number of dimensions for data
                    iniparams[3] = len(j)
                    iniparams[4] = filelength
                    iniparams.tofile(fidout)

                    fidin = [open(k, 'r') for k in j]

                    c = 1
                    br = 0
                    while True:
                        towrite = []
                        for k in fidin:
                            towrite.append(np.fromfile(k, 'int16', nsamples))
                            if towrite[0].size == 0:
                                br = 1
                                break
                        if br:
                            break
                        towrite = np.asarray(towrite, 'int16')
                        towrite = towrite.flatten('F')
                        towrite.tofile(fidout)
                        print('{} samples processed...'.format(c * nsamples))
                        c += 1


def batch_save_mountainsort_results(outfolder, folder=r'C:\Users\Jermyn\BansheeShared\MSResults', savemethod='scipymat',
                                    probefolder=r'C:\Users\Jermyn\BansheeShared\MSDefaultFiles', triggerfolder=(),
                                    combineopt=False, order=()):

    subfolders = glob.glob(os.path.join(folder, '*-firings'))
    basefolders = [os.path.basename(i) for i in subfolders]

    if combineopt:
        folderlist = match_multiprobe_MS_folders(basefolders, order)
    else:
        folderlist = [[i] for i in basefolders]


    for i in folderlist:

        intdepth = []
        spk = []
        allprobe = ''
        alldepth = ''

        for j in i:

            basefile = re.search('^\S+(?=-firings)', j).group(0)

            if not os.listdir(os.path.join(folder, j)):
                print('Output files missing from {}. Skipping...'.format(j))
                continue

            else:

                print('Processing {}...'.format(j))

                depth = re.search('(?<=\d-)\w+(?=um-)', basefile).group(0)
                # parse depth if necessary
                if not depth.isdigit():
                    temp = re.findall('\d{2,5}', depth)
                    tempdepth = [int(i) for i in temp]
                    if len(temp) == 1:
                        temp = tempdepth[0]
                    intdepth.append(temp)
                else:
                    intdepth = int(depth)

                mdafile = os.path.join(folder, j, 'firings_uncurated.mda')
                jsonfiles = os.path.join(folder, j, 'cluster_metrics.json')
                # pairjsonfile = os.path.join(folder, j, 'pair_metrics.json')
                # try:
                probe = re.search('(?<=min-)\S+(?=-fs)', basefile).group(0)
                # except AttributeError:
                #     probe = re.search('(?<=combined-)\S+(?=-fs)', basefile).group(0)

                probefile = probe + '.csv'
                csvfile = os.path.join(probefolder, '{}'.format(probefile))


                probechanstr = re.findall('\d{1,2}x\d{1,3}', probefile)
                chancount = int(np.prod(np.asarray([re.split('x',i) for i in probechanstr], dtype=int), axis=1))

                temp = process_mountainsort_results(mdafile, jsonfiles, csvfile, intdepth, chancount)

                if combineopt:
                    spk.append(temp)
                    if not allprobe:
                        allprobe = allprobe + probe
                    else:
                        allprobe = allprobe + '-' + probe
                    alldepth += depth

                else:
                    spk = temp
                    allprobe = probe
                    alldepth = depth

        if combineopt:
            part1 = re.search('^\S+site\d{1,2}-', basefile).group(0)
            part2 = re.search('-\d{1,3}db\S+min-', basefile).group(0)
            part3 = re.search('-fs\d{4,6}', basefile).group(0)
            alldepth += 'um'
            basefile = part1 + alldepth + part2 + allprobe + part3
            spk = [j for i in spk for j in i]

        if savemethod == 'json':
            outfile = os.path.join(outfolder, basefile + '-pydict.json')
        elif savemethod == 'scipymat':
            outfile = os.path.join(outfolder, basefile + '-pydict.mat')
        else:
            raise ValueError("Use either json or scipymat as the save method!")

        dictsave = {
            'exp': re.search('^\d{6}_\d{6}(?=-site)', basefile).group(0),
            'site': float(re.search('(?<=-site)\d{1,2}(?=-)', basefile).group(0)),
            'stim': re.search('(?<=db-)\w+(?=-)', basefile).group(0),
            'depth': float(intdepth),
            'atten': re.search('(?<=um-)\w+db', basefile).group(0),
            'fs': float(re.search('(?<=-fs)\d{4,6}$', basefile).group(0)),
            'probe': allprobe,
            'stimlength': re.search('(?<=-)\d{1,3}min(?=-)', basefile).group(0),
            'spk': spk
        }


        if savemethod == 'json':

            with open(outfile, 'w') as fidout:
                json.dump(dictsave, fidout, indent=4)

        elif savemethod == 'scipymat':
            sio.savemat(outfile, {'spk':dictsave})

        if triggerfolder:

            trigname = basefile + r'-ADC-00.mat'
            trigfile = []

            for paths, _, files in os.walk(triggerfolder):
                if trigname in files:
                    trigfile.append(os.path.join(paths, trigname))

            if len(trigfile) == 1:
                shutil.copy(trigfile[0], os.path.join(outfolder, trigname))
            elif len(trigfile) == 0:
                warnings.warn('{} has no trigger file!'.format(basefile))
            else:
                raise Exception('More than one trigger file found!')

    # if combineopt:



#
# def get_recording_params(subfolder, inputfolder=):
#
#     foldername = os.path.basename(subfolder)
#
#
#     inputfolder = r'../MSInput'
#     inputmdafile = glob.glob(os.path.join(inputfolder, exp + '-' + site + '*-' + stim + '-*'))
#
#     assert len(inputmdafile) == 1
#
#     inputmdafile = os.path.basename(inputmdafile[0])
#
#     return re.search('^\S+(?=-mountainsort.mda)', inputmdafile).group(0)



def process_mountainsort_results(mdafile, jsonfiles, csvfile, depth, chancount):

    spkmat = readmda(mdafile)
    spkdata = parse_mda_results(spkmat, csvfile, depth, chancount)

    paramtemp = []
    if isinstance(jsonfiles, str):
        jsonfiles = [jsonfiles]


    for jf in jsonfiles:
        with open(jf, 'r') as jfile:
            temp = json.load(jfile)
            paramtemp.append(temp['clusters'])

    params = []

    for pt in paramtemp:
        if not params:
            params = pt
        else:
            [i['metrics'].update(j['metrics']) for i,j in zip(params,pt)]

    # to prevent bugs in saving to .mat files (None cannot be saved)

    for i in params:

        if not i['metrics']['peak_noise']:
            i['metrics']['peak_noise'] = 0
            i['metrics']['peak_snr'] = 0

    if len(params) == len(spkdata):
        spk = [{**spkdata[i], **params[i]['metrics']} for i in range(len(params))]
    else:
        idx = [i['unit'] - 1 for i in spkdata]
        spk = [{**spkdata[i], **params[j]['metrics']} for i, j in zip(range(len(spkdata)), idx)]

    try:
        spk = sort_mountainsort_data_by_position(spk)
    except:
        spk = []

    return spk

def sort_mountainsort_data_by_position(spk):

    pos = [i['position'] for i in spk]
    xpos = np.asarray([i[0] for i in pos])
    ypos = np.asarray([i[1] for i in pos])
    yunique = np.unique(ypos)
    sortidx = np.argsort(ypos)
    ysorted = ypos[sortidx]
    xsorted = xpos[sortidx]

    fsortidx = []
    for i in yunique:
        yidx = ysorted == i
        xidx = np.argsort(xsorted[yidx])
        tempidx = sortidx[yidx]
        fsortidx.append(tempidx[xidx])

    fsortidx = np.concatenate(fsortidx).tolist()

    return [spk[i] for i in fsortidx]




def readmda(mdafile):

    with open(mdafile, 'r') as fid:
        mdaparams = np.fromfile(fid, 'int32', 3)
        code = mdaparams[0]
        ndim = mdaparams[2]

        S = np.empty([ndim])

        for i in range(ndim):
            S[i] = np.fromfile(fid, 'int32', 1)

        S = tuple(map(int, S))

        N = int(np.prod(S))

        if code == -4:
            A = np.fromfile(fid, 'int16', N)
        elif code == -5:
            A = np.fromfile(fid, 'int32', N)
        elif code == -7:
            A = np.fromfile(fid, 'd', N)

        return np.reshape(A, S, 'F')


def parse_mda_results(resultsmat, csvfile, depth, chancount):

    resultsmat[0] = resultsmat[0] - 1 #change channel number to start from 0, which is intan's format

    spkdata = []
    units = np.unique(resultsmat[2])

    position = get_position_from_geom(csvfile, depth, chancount)

    for i in units:
        idx = np.where(resultsmat[2,:] == i)
        temp = {
            'unit': int(i),
            'chan': int(resultsmat[0, idx[0][0]]),
            'spktimes_samp': tuple(map(int,tuple(resultsmat[1, idx][0]))),
            # 'amplitude': tuple(resultsmat[3, idx][0]),
            'position': position[int(resultsmat[0, idx[0][0]])]
        }
        spkdata.append(temp)

    return spkdata


def get_position_from_geom(csvfile, depth, chancount):

    with open(csvfile, 'r') as cfile:
        position = list(csv.reader(cfile))

    if isinstance(depth, int):
        assert(isinstance(chancount, int))
        assert(len(position) == chancount)
        depthvec = np.tile(depth, chancount)
    elif isinstance(depth, list):
        assert(len(depth) == len(chancount))
        assert(len(position) == np.sum(chancount))
        depthvec = np.asarray([np.tile(i, j) for i,j in zip(depth, chancount)]).flatten()
    else:
        raise Exception('Depth must be an integer or a list of integers.')

    position = [list(map(int, i)) for i in position]
    position = [[i[0], j - i[1]] for i,j in zip(position, depthvec)]

    return position

def match_multiprobe_MS_folders(folders, order=()):

    exp = [re.search('^\d{6}_\d{6}', i).group(0) for i in folders]
    stim = [re.search('(?<=db-)\S+min', i).group(0) for i in folders]

    expidx = [np.where(np.asarray(exp) == x) for x in np.unique(exp)]
    stimidx = [np.where(np.asarray(stim) == x) for x in np.unique(stim)]

    matchedidx = [np.intersect1d(x,y) for x in expidx for y in stimidx]
    folderarray = []

    for i in range(len(matchedidx)):
        temp = []
        for j in range(matchedidx[i].size):
            temp.append(folders[matchedidx[i][j]])
        folderarray.append(temp)

    if order:

        for i in range(len(folderarray)):
            id = [re.search('[a-zA-Z]+(?=\d{2,6}um)', foldername).group(0) for foldername in folderarray[i]]
            idx = [id.index(o) for o in order]
            folderarray[i] = [folderarray[i][ii] for ii in idx]


    return folderarray


#
# def raw2spykingcircus(folders, stim, prbfile, batfileopt=1):
#
#     """
#     # folder: list of folders (normally separated by different sites)
#     # stim: ['rn1', 'rn4', etc.]
#     """
#
#     #check if prbfile exists
#     assert os.path.isfile(r'I:\Data\SpykingCircusMasterFiles\\{}'.format(prbfile)), '.prb file does not exist'
#
#     if batfileopt == 1:
#         batfilefolders = []
#
#     for folds in folders:
#
#         os.chdir(folds)
#         nsamples = 1000000
#
#         for i in stim:
#
#             files = glob.glob('*-site*-*um-*db-{}-fs*-*-*.raw'.format(i))
#
#             #remove trig file
#             trigidx = next((i for i, j in enumerate([re.findall('ADC', k) for k in files]) if j), len(files))
#             del files[trigidx]
#
#             exp = re.search('^\d{6}_\d{6}', files[0]).group(0)
#
#             #get subfolder name
#             subfoldername = re.sub('(?<=db_)\w+(?=_\d{6}_\d{6})', i, os.getcwd().split('\\')[-1])
#             sfexp = re.search('\d{6}_\d{6}$', subfoldername).group(0)
#
#             #correct subfolder name timestamp if necessary
#             if sfexp != exp:
#                 subfoldername = re.sub(sfexp, exp, subfoldername)
#
#             #create subfolder if it doesn't exist
#             if not os.path.isdir('.\\{}'.format(subfoldername)):
#                 os.mkdir('.\\{}'.format(subfoldername))
#
#             #get outfile name
#             filebase = re.match('^\S+(?=-[ABCD]-\d{3}.raw)', files[0])
#             outfile = ('\\{}\\{}_spyking_circus.int16'.format(subfoldername, filebase.group(0)))
#
#             #skip existing files
#             if os.path.isfile(outfile):
#                 print('{} already exists! Skipping...'.format(outfile))
#                 continue
#
#             if batfileopt == 1:
#                 batfilefolders.append(os.getcwd() + outfile)
#
#             outfile = '.' + outfile
#
#             # if not os.path.exists('./spykingcircus'):
#             #     os.mkdir('./spykingcircus/')
#
#             fidin = [open(j, 'r') for j in files]
#             fidout = open(outfile, 'w')
#
#             c = 1
#             br = 0
#             while True:
#                 towrite = []
#                 for j in fidin:
#                     towrite.append(np.fromfile(j, 'int16', nsamples))
#                     if towrite[0].size == 0:
#                         br = 1
#                         break
#                 if br:
#                     break
#                 towrite = np.asarray(towrite, 'int16')
#                 towrite = towrite.flatten('F')
#                 towrite.tofile(fidout)
#                 print('{} samples processed...'.format(c * nsamples))
#                 c += 1
#
#             # shutil.copy(r'I:\Data\SpykingCircusMasterFiles\config.params', '.\\{}\\{}_spyking_circus.params'.format
#             # (subfoldername, filebase.group(0)))
#             subfolder = '.\\{}'.format(subfoldername)
#             paramfile = shutil.copy(r'I:\Data\SpykingCircusMasterFiles\config.params', subfolder)
#
#             # copy paramfile
#             with open(paramfile) as f:
#                 paramtext = f.readlines()
#
#             # find and add probefile
#             mappingidx = next((i for i, j in enumerate(['mapping' in k for k in paramtext]) if j), len(paramtext))
#
#             paramtext[mappingidx] = re.sub('(?<==)\s', r' I:\Data\SpykingCircusMasterFiles\\{}'.format(prbfile),
#                                            paramtext[mappingidx])
#
#             # write paramfile with correct name
#             outparamfile = '.\\{}\\{}_spyking_circus.params'.format(subfoldername, filebase.group(0))
#
#             with open(outparamfile, 'w') as file:
#                 file.writelines(paramtext)
#
#             # delete default paramfile
#             os.remove(paramfile)
#
#
#     if batfileopt == 1:
#         os.chdir(r'I:\Python\Python Codes')
#         create_spyking_circus_bat_file(batfilefolders)
#
#     return
#
#
# def create_spyking_circus_bat_file(int16files):
#     with open('spyking_circus_batch_run.bat', 'w') as batfile:
#         batfile.write('@echo off\n')
#         for i in int16files:
#             batfile.write('spyking-circus ' + i + '\n')
#     return