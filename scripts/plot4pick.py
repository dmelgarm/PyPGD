def plot4pick(k,event,ID):
    
    from matplotlib import pyplot as plt
    from obspy import read
    from numpy import mean
    
    plt.close('all')
    sta=event.filtered_stations[k]
    st=event.waveforms[k]
    
    plt.figure(figsize=(14,3))
    plt.plot(st.times(),st.data-mean(st.data),zorder=0)
    plt.scatter(event.picks[k],0,c='r');
    plt.title(sta)
    
    plt.show()
    
    
def write4picking(event,folder_out):
    
    for k in range(len(event.filtered_stations)):
        tr=event.waveforms[k]
        out=folder_out+tr.stats.station+'.'+tr.stats.channel+'.sac'
        tr.write(out,format='SAC')
        
        
def write_picks(folder):
    
    from obspy import read
    from glob import glob
    from numpy import nan
    
    files=glob(folder+'*HNZ*')
    
    for k in range(len(files)):
        st=read(files[k])
        station=st[0].stats.station
        try:
           pick=st[0].stats['sac']['a']
        except:
           pick=nan 
    
        print station,pick
    