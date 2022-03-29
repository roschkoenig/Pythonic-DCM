import mne
import numpy as np
import networkx as nx
import community
import matplotlib.pyplot as plt
import pyCM
import scipy
import os

#=================================================================================
def dcm_spec(dcm_path, labels, time, datatype = 'CSD', model = 'CMC', spatial = 'LFP',
             freq_range = [1,30], Fwd = None, Bwd = None, Slf = None):
# Input: key parameters for DCM inversion
# Output: specified DCM structure as dictionary
#=================================================================================
    # Settings for Dynamic Causal Modelling
    #==========================================================================
    DCM = {'options':{}, 'A':[], 'B':[], 'C':[], 'M':{'dipfit':{}}}
    DCM['name'] = dcm_path
             
    # DCM modelling settings
    #--------------------------------------------------------------------------
    DCM['options']['analysis'] = datatype
    DCM['options']['model']    = model
    DCM['options']['spatial']  = spatial
    DCM['Sname']               = labels
             
    # Signal and preprocessing settings
    #--------------------------------------------------------------------------
    DCM['options']['Tdcm']     = [time[0], time[-1]*1000]   # start and end time in miliseconds
    DCM['options']['Fdcm']     = freq_range  # frequency range 
    DCM['options']['D']        = 1       # frequency bin, 1 = no downsampling
    DCM['options']['Nmodes']   = 8       # cosinde reduction components used 
    DCM['options']['han']      = 0       # no Hanning
    DCM['options']['trials']   = 1       # ! which condition will be modelled; indexing starting with 1
    
    # Define possible connections
    #--------------------------------------------------------------------------
    if not Fwd: Fwd  = np.triu(np.ones(len(DCM['Sname'])),1) # upper triangle = all possible forward connections 
    if not Bwd: Bwd  = np.tril(np.ones(len(DCM['Sname'])),-1) # lower triangle = all possible backward connections 
    if not Slf: Slf  = np.diag(np.ones(len(DCM['Sname'])))   # recurrent self connections along diagonal
    
    # Define model architecture (A), conditional effects (B), and input (C)
    #--------------------------------------------------------------------------
    DCM['A'].append(Fwd+Bwd)  # equivalent to A{1} in Matlab: forward connections
    DCM['A'].append(Fwd+Bwd)  # equivalent to A{2} in Matlab: backward connections
    DCM['A'].append(Slf)    # equivalent to A{3} in Matlab: self connections

    DCM['B'] = []         # no conditional effects on connectivity 

    DCM['C'] = []         # no specified input
    
    # Package up options so that data can be estimated independently
    #---------------------------------------------------------------------------
    DCM['M']['dipfit']['model'] = DCM['options']['model']
    DCM['M']['dipfit']['spatial'] = DCM['options']['spatial']
    
    return(DCM)

#=================================================================================
def dataload(filename, elect=[], chans=[], seg_type = 'continuous', window_secs = 0,
             window_number=1,downsample=250):
# Input: filename from BIDS datastructure
# Output: epoched SEEG object epochs*channels*samples
#=================================================================================

    # identify good channels from whole electrode list
    #------------------------------------------------------------------------------
    goodels       = elect[elect['cortical'] == True]
    badrows       = [rowid for rowid,row in chans.iterrows() 
                     if goodels[goodels['name'] == row['name']].empty]
    goodrows      = [rowid for rowid,row in chans.iterrows() 
                     if not goodels[goodels['name'] == row['name']].empty]

    # Load raw data (memory map, not preload) 
    #------------------------------------------------------------------------------
    raw              = mne.io.read_raw_edf(filename)
    raw.info['bads'] = list(chans['name'][badrows])
    Nsampls          = raw.n_times
    Fs               = raw.info['sfreq']

    # Sample from continuous raw data according to sampling scheme
    #------------------------------------------------------------------------------
    if seg_type == 'uniform':
        winsamps = window_secs * Fs
        events   = np.multiply([range(window_number)],int((Nsampls - winsamps)/window_number))
        eventarr = np.zeros([events.shape[1],3])
        eventarr[:,0] = events
        eventarr[:,2] = np.ones(events.shape)

    if seg_type == 'continuous':
        eventarray = np.array([0,0,1])

    epochs = mne.Epochs(raw, eventarr.astype('int'), tmin=0, tmax=window_secs,baseline=None,
                        preload=True)
    if downsample: 
        epochs.resample(downsample)
        ds_statement = 'Data resampled to '+str(downsample)+'Hz'
    else: 
        ds_statement = 'Data maintained at original sampling frequency'
    dat = epochs.get_data()[:,goodrows,:]
    lbl = [epochs.ch_names[gr] for gr in goodrows]
    tim = epochs.times
    Fs  = Fs

    # Pack up for export
    #------------------------------------------------------------------------------
    D = {'data':dat, 'labels':lbl, 'time':tim, 'Fs':Fs}
    
    # Status update
    #------------------------------------------------------------------------------
    print('------------------\nData extracted \nSampling: ' +seg_type+ '\n'+
          'Epochs: ' +str(window_number)+ ' (' +str(window_secs)+' seconds each)\n'+
          ds_statement)
    return(D)

#=================================================================================
def node_selection(seeg, method='louvain', output='single-node', cutoff=10, reps=50):
# Input: epoched SEEG dataset extracted through dataload epochs*channels*samples
# Output: reduced dataformat with <10 channels epochs*reduced-channels*samples
#=================================================================================
    Nreps = seeg['data'].shape[0]
    Nchans = seeg['data'].shape[1]
    
    #----------------------------------------------------------------------------
    if method == 'louvain': # <- Performs partitioning
    # 
    # This method applies community detection the the covariance matrix 
    # It will identify a consensus partition based on repeat partitonig
    #----------------------------------------------------------------------------
        cov = np.zeros((Nreps,Nchans,Nchans))
        for di in range(Nreps):
            d   = seeg['data'][di,:,:]
            cov[di,:,:] = np.cov(d)

        mcov  = np.mean(cov,axis=0)
        psize = 0
        G     = nx.from_numpy_array(mcov)
        for resolution in range(1, 10, 1):
            partition = community.best_partition(G,resolution=resolution/10)
            this_psize = len(set(partition.values()))
            if (this_psize > psize) & (this_psize <= cutoff):
                psize = this_psize
                final_resolution = resolution/10
                final_partition = partition
        print('Resolution fixed at '+str(final_resolution))

        # Repeat Partitioning for robustness
        #--------------------------------------------------------------------------
        for i in range(reps):
            Acomms    = np.zeros((reps,mcov.shape[0],mcov.shape[1]))
            partition = community.best_partition(G,resolution=final_resolution)
            p         = [partition[pi] for pi in range(len(partition))]
            for pi in range(len(p)):
                samesies = np.where(np.array(p) == p[pi])[0]
                Acomms[i,pi,samesies] = 1

        cG = nx.from_numpy_array(np.mean(Acomms,axis=0))
        partition = community.best_partition(cG,resolution=final_resolution)
        print('Consensus partition with '+str(len(set(partition.values())))+' communities')
        A = cov
        
    #----------------------------------------------------------------------------
    if output=='single-node': # <- Selects single nodes from partition
    # 
    # This method will identify the node (per EEG segment) that correlates most 
    # with the mean trace and select this as 'representative node' for further
    # analysis 
    #----------------------------------------------------------------------------
        pvec = [partition[pi] for pi in range(len(partition))]
        ps = np.unique(pvec)
        nodelist = []
        for p in ps:
            nodes = np.where(pvec == p)[0]
            if len(nodes) == 1: node = nodes[0]
            else: 
                dat  = seeg['data'][:,nodes,:]
                mdat = np.mean(dat,axis=1) 
                mcc  = []
                for di in range(len(nodes)):
                    cc=[]
                    for si in range(seeg['data'].shape[0]):
                        cc.append(np.corrcoef(dat[si,di,:],mdat[si])[0,1])
                    mcc.append(np.mean(cc))
                node = nodes[np.argmax(mcc)]
            nodelist.append(node)

    return(A, partition,nodelist)

#=================================================================================
def plot_by_partition(A,partition,ax=None):
# Input: n*n adjacency matrix and n-dimensional (dictionary) partition 
# Output: 2d matrix plot of adjacency matrix ordered by partition
#=================================================================================
    pvec    = [partition[pi] for pi in range(len(partition))]
    sorting = np.argsort(pvec)
    srtd    = np.array(pvec)[sorting]

    # Identify transitions between partitions
    #----------------------------------------------------------------------------
    curval    = srtd[0]
    stp       = -1
    segments  = []
    stepcount = 0

    for p in srtd:
        if curval != p:
            srt = stp+1
            stp = stepcount-1
            segments.append((srt,stp))
        curval = p
        stepcount = stepcount + 1
    segments.append((stp,stepcount-1))
    
    # Sort the adjacency matrix by partition
    #----------------------------------------------------------------------------
    tA = A[sorting,:]
    tA = tA[:,sorting]
    
    # Plotting
    #----------------------------------------------------------------------------
    if ax==None: fig,ax = plt.subplots(1,1)
    
    # Plot adjacency matrix
    image = ax.imshow(tA,cmap='gray')
    ax.set_xlabel('nodes'), ax.set_xticks([]), ax.set_yticks([])
    ax.set_ylabel('nodes'), ax.set_xticks([]), ax.set_yticks([])
    plt.colorbar(image, ax=ax, shrink=0.7)
    
    # Plot coloured lines to indicate segments
    colors = plt.get_cmap(name='tab10', lut=len(segments))
    
    i = 0
    for s in segments:
        ax.plot([0,0],[s[0],s[1]], color=colors(i), linewidth=10)
        ax.plot([s[0],s[1]],[0,0], color=colors(i), linewidth=10)
        i = i+1

def start_matlab(F):
    import matlab.engine
    M = matlab.engine.start_matlab()
    M.addpath(F['code']+os.sep+'matlab')
    print('Matlab Engine started')
    return M
        
# List of outstanding issues
# 
# - Louvain community detection not particularly robust to repetition
# - for node selection, a virtual electrode approach could be performed 
# - loading from .mat files adds weird empty dimensions to arrays