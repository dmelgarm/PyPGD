


def knet2stream(knetfile):
    '''
    convert KNET ascii format to sac
    '''
    
    from obspy import Stream,Trace,UTCDateTime
    from numpy import genfromtxt,size
    
    #Parse headers
    f=open(knetfile,'r')
    
    while True:
        line=f.readline()
        if 'Record Time' in line:
            tstart=line.split()[2]+'T'+line.split()[3]
            tstart=tstart.replace('/','-')
            tstart=UTCDateTime(tstart)
        if 'Station Code' in line:
            station=line.split()[-1]
        if 'Sampling Freq' in line:
            fs=float(line.split()[-1].replace('Hz',''))
            dt=1./fs
        if 'Dir.' in line:
            component=line.split()[-1]
            if component=='E-W':
                channel='HNE'
            if component=='N-S':
                channel='HNN'
            if component=='U-D':
                channel='HNZ'
        if 'Scale Factor' in line:
            scale1=line.split()[-1].split('/')[0].replace('(gal)','')
            scale2=line.split()[-1].split('/')[1]
            scale=float(scale1)/float(scale2)
            break
    f.close()
    
    #Now read data
    data=genfromtxt(knetfile,skip_header=17,skip_footer=1)
    data=data.reshape(1,size(data)).squeeze()
    data=data*scale
    
    #Knet is cgs
    data=data/100.
    
    #Output
    st=Stream(Trace())
    st[0].data=data
    st[0].stats.delta=dt
    st[0].stats.starttime=tstart
    st[0].stats.station=station
    st[0].stats.channel=channel
    
    return st
    
    
    
def process1folder(infolder,outfolder):
    '''
    Convert all data from a fodler
    '''
    
    from glob import glob
    
    #Which files?
    
    files=glob(infolder+'*EW')
    files+=glob(infolder+'*NS')
    files+=glob(infolder+'*UD')
    
    for k in range(len(files)):
        print files[k]
        st=knet2stream(files[k])
        outname=outfolder+st[0].stats.station+'.'+st[0].stats.channel+'.sac'
        st.write(outname,format='SAC')