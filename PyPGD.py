    
class Event(object):
    '''
    Things related to a particular earthquakes
    '''

    def __init__(self,path_to_data,path_to_station_file,hypocenter,hypo_time,duration):
        
        '''
        Initalize the class
        '''        
        
        from numpy import genfromtxt
        
        self.hypocenter=hypocenter
        self.hypo_time=hypo_time
        self.path_to_data=path_to_data+'/proc/'
        self.path_to_station_file=path_to_station_file
        self.duration=duration
        self.station_locations=genfromtxt(self.path_to_station_file,usecols=[1,2])
        self.station_names=genfromtxt(self.path_to_station_file,usecols=0,dtype='S')
        self.Nsta=len(self.station_names)
        
    
    def get_travel_times(self):
        '''
        Get direct P and S to each station
        '''
        
        from obspy.taup import TauPyModel
        from numpy import float64,ones
        from obspy.geodetics import locations2degrees
        from datetime import timedelta
        

        self.ptime=1e6*ones(self.Nsta)
        self.stime=1e6*ones(self.Nsta)        
                        
        for k in range(self.Nsta):
        
            #tauP class
            velmod = TauPyModel()
        
            #Get station to hypocenter delta distance
            delta=locations2degrees(self.station_locations[k,1],self.station_locations[k,0],self.hypocenter[1],self.hypocenter[0])
            
            
            
            #Get travel times
            arrivals = velmod.get_travel_times(source_depth_in_km=self.hypocenter[2],
                        distance_in_degree=delta,phase_list=['P','Pn','S','Sn','p','s'])
            
            #Determine P and S arrivals
            for kphase in range(len(arrivals)):
                if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
                    if arrivals[kphase].time<self.ptime[k]:
                        self.ptime[k]=arrivals[kphase].time
                if 'S' == arrivals[kphase].name or 's' == arrivals[kphase].name or 'Sn' == arrivals[kphase].name:
                    if arrivals[kphase].time<self.stime[k]:
                        self.stime[k]=arrivals[kphase].time
                        
    
    def get_SNR(self,t_prior=60,source_durations=1.5):
        '''
        Get SNR considering t_prior seconds before p-wave and a
        total waveofrm of source_durations length
        
        '''
        
        from numpy import zeros,log10
        from obspy import read
        from datetime import timedelta
        from obspy import UTCDateTime
        from string import rjust
        
        self.SNR_north=zeros(self.Nsta)  
        self.SNR_east=zeros(self.Nsta)   
        self.SNR_up=zeros(self.Nsta) 
        
        hypo_time=UTCDateTime(self.hypo_time)
        
        for k in range(self.Nsta):
            #Are the waveform slices empty?
            empty=False
            
            path=self.path_to_data
            
            ptime=self.ptime[k]

            tp=hypo_time+timedelta(seconds=ptime)  
            t0=tp-timedelta(seconds=t_prior)
            tfinal=tp+timedelta(seconds=source_durations*self.duration)          
                                      
            #read waveforms
            try:
                sta=rjust(str(self.station_names[k]),4,'0')
                n=read(path+sta+'.LXN.sac')
                e=read(path+sta+'.LXE.sac')
                z=read(path+sta+'.LXZ.sac')
            except:
                print 'ERROR: '+self.station_names[k]+'.LXN.sac not found'
                return
            
            #Get pre-event portions (from t0 to tp)
            n_pre=n.copy()
            e_pre=e.copy()
            z_pre=z.copy()
            
            n_pre[0].trim(starttime=t0,endtime=tp)
            e_pre[0].trim(starttime=t0,endtime=tp)
            z_pre[0].trim(starttime=t0,endtime=tp)
            
            if n_pre[0].stats.npts==0 or e_pre[0].stats.npts==0 or z_pre[0].stats.npts==0:
                empty=True
            else:
                mean_n=n_pre[0].data.mean()
                mean_e=e_pre[0].data.mean()
                mean_z=z_pre[0].data.mean()
                
                n_pre[0].data=n_pre[0].data-mean_n
                e_pre[0].data=e_pre[0].data-mean_e
                z_pre[0].data=z_pre[0].data-mean_z
                
                #get power
                Pn_pre=sum(n_pre[0].data**2)/n_pre[0].stats.npts
                Pe_pre=sum(e_pre[0].data**2)/e_pre[0].stats.npts
                Pz_pre=sum(z_pre[0].data**2)/z_pre[0].stats.npts
            
            #Trim post P-arrival signal
            n[0].trim(starttime=tp,endtime=tfinal)
            e[0].trim(starttime=tp,endtime=tfinal)
            z[0].trim(starttime=tp,endtime=tfinal)
            
            
            if n[0].stats.npts==0 or e[0].stats.npts==0 or z[0].stats.npts==0:
                empty=True
            else:
                n[0].data=n[0].data-mean_n
                e[0].data=e[0].data-mean_e
                z[0].data=z[0].data-mean_z
                
                #get power
                Pn=sum(n[0].data**2)/n[0].stats.npts
                Pe=sum(e[0].data**2)/e[0].stats.npts
                Pz=sum(z[0].data**2)/z[0].stats.npts
            
            #And finaly SNR
            if empty==False:
                self.SNR_north[k]=Pn/Pn_pre
                self.SNR_east[k]=Pe/Pe_pre
                self.SNR_up[k]=Pz/Pz_pre
            
        
    def get_hypocentral_distance(self):
        '''
        Compute distance from stations to hypocenter, considers depth
        '''
        
        from numpy import zeros
        from obspy.geodetics import gps2dist_azimuth
        
        self.hypo_distance_in_km=zeros(self.Nsta)
        
        for k in range(self.Nsta):
            lat1=self.station_locations[k,1]
            lon1=self.station_locations[k,0]
            lat2=self.hypocenter[1]
            lon2=self.hypocenter[0]
            d,Az,Baz=gps2dist_azimuth(lat1, lon1, lat2, lon2)
            d=((d/1000.)**2+self.hypocenter[2]**2)**0.5
            self.hypo_distance_in_km[k]=d
 
               
    def get_PGD(self,t_prior=60,source_durations=2.0):
        '''
        Get SNR considering t_prior seconds before p-wave and a
        total waveofrm of source_durations length
        
        '''
        
        from numpy import zeros,log10
        from obspy import read
        from datetime import timedelta
        from obspy import UTCDateTime
        from string import rjust
        
        self.PGD_north=zeros(self.Nsta)  
        self.PGD_east=zeros(self.Nsta)   
        self.PGD_up=zeros(self.Nsta) 
        self.PGD_cart=zeros(self.Nsta) 
        
        hypo_time=UTCDateTime(self.hypo_time)
        
        for k in range(self.Nsta):
            
            #Are the waveform slices empty?
            empty=False
            
            path=self.path_to_data
            
            ptime=self.ptime[k]

            tp=hypo_time+timedelta(seconds=ptime)  
            t0=tp-timedelta(seconds=t_prior)
            tfinal=tp+timedelta(seconds=source_durations*self.duration)          
                                      
            #read waveforms
            try:
                sta=rjust(str(self.station_names[k]),4,'0')
                n=read(path+sta+'.LXN.sac')
                e=read(path+sta+'.LXE.sac')
                z=read(path+sta+'.LXZ.sac')
            except:
                print 'ERROR: '+self.station_names[k]+'.LXN.sac not found'
                return
            
            #Get pre-event portions (from t0 to tp)
            n_pre=n.copy()
            e_pre=e.copy()
            z_pre=z.copy()
            
            n_pre[0].trim(starttime=t0,endtime=tp)
            e_pre[0].trim(starttime=t0,endtime=tp)
            z_pre[0].trim(starttime=t0,endtime=tp)
            
            if n_pre[0].stats.npts==0 or e_pre[0].stats.npts==0 or z_pre[0].stats.npts==0:
                empty=True
            else:
                #Get pre-event mean to remove from waveform
                mean_n=n_pre[0].data.mean()
                mean_e=e_pre[0].data.mean()
                mean_z=z_pre[0].data.mean()
                
            
            #Trim post P-arrival signal
            n[0].trim(starttime=tp,endtime=tfinal)
            e[0].trim(starttime=tp,endtime=tfinal)
            z[0].trim(starttime=tp,endtime=tfinal)
            
            if n[0].stats.npts==0 or e[0].stats.npts==0 or z[0].stats.npts==0:
                empty=True
            else:
                n[0].data=n[0].data-mean_n
                e[0].data=e[0].data-mean_e
                z[0].data=z[0].data-mean_z
                
                #Construct "cartesian" sum waveform
                c=n.copy()
                c[0].data=(n[0].data**2+e[0].data**2+z[0].data**2)**0.5
                
            if empty==False: #there was enough data
                #get PGD
                self.PGD_north[k]=abs(n[0].data).max()
                self.PGD_east[k]=abs(e[0].data).max()
                self.PGD_up[k]=abs(z[0].data).max()
                self.PGD_cart[k]=abs(c[0].data).max()
            
                           
    def get_PGD_with_time(self,tinit='ptime',tfinal=200):
        '''
        Get PGD for every station as a function of time
        
        '''
        
        from numpy import where,c_,array,nan
        from obspy import read
        from obspy import UTCDateTime
        from string import rjust
        from scipy import maximum
        
        self.PGD_north_with_time=[]
        self.PGD_east_with_time=[]   
        self.PGD_up_with_time=[] 
        self.PGD_cart_with_time=[]
        
        
        hypo_time=UTCDateTime(self.hypo_time)
        
        for k in range(self.Nsta):
            
            #Are the waveform slices empty?
            empty=False
            
            path=self.path_to_data
            
            if tinit=='ptime':
                tstart=self.ptime[k]
            else:
                tstart=self.stime[k]        
                                      
            #read waveforms
            try:
                sta=rjust(str(self.station_names[k]),4,'0')
                n=read(path+sta+'.LXN.sac')
                e=read(path+sta+'.LXE.sac')
                z=read(path+sta+'.LXZ.sac')
            except:
                print 'ERROR: '+self.station_names[k]+'.LXN.sac not found'
                return
            
            #Trim from event origin to end
            n[0].trim(starttime=hypo_time,endtime=hypo_time+tfinal)
            e[0].trim(starttime=hypo_time,endtime=hypo_time+tfinal)
            z[0].trim(starttime=hypo_time,endtime=hypo_time+tfinal)
            
            #FIND TIMES BEFORE TSTART (p or s time)
            inorth=where(n[0].times()<tstart)[0]
            ieast=where(e[0].times()<tstart)[0]
            iup=where(z[0].times()<tstart)[0]
            
            if n[0].stats.npts==0 or e[0].stats.npts==0 or z[0].stats.npts==0:
                empty=True
            else:
                #Get pre-event mean to remove from waveform
                mean_n=n[0].data[inorth].mean()
                mean_e=e[0].data[ieast].mean()
                mean_z=z[0].data[iup].mean()
                
            
            if n[0].stats.npts==0 or e[0].stats.npts==0 or z[0].stats.npts==0:
                empty=True
            else: #de-mean and make things before arrival =0 (no PGD possible)
                n[0].data=n[0].data-mean_n
                e[0].data=e[0].data-mean_e
                z[0].data=z[0].data-mean_z
                n[0].data[inorth]=0
                e[0].data[ieast]=0
                z[0].data[iup]=0
                
                #Construct "cartesian" sum waveform
                c=n.copy()
                c[0].data=(n[0].data**2+e[0].data**2+z[0].data**2)**0.5
                
            if empty==False: #there was enough data
                #get PGD include time vector
                north=c_[n[0].times(),maximum.accumulate(abs(n[0].data))]
                east=c_[e[0].times(),maximum.accumulate(abs(e[0].data))]
                up=c_[z[0].times(),maximum.accumulate(abs(z[0].data))]
                cart=c_[c[0].times(),maximum.accumulate(abs(c[0].data))]
                
                self.PGD_north_with_time.append(north)
                self.PGD_east_with_time.append(east)
                self.PGD_up_with_time.append(up)
                self.PGD_cart_with_time.append(cart) 
                
            else:                                
                self.PGD_north_with_time.append(array([nan]))
                self.PGD_east_with_time.append(array([nan]))
                self.PGD_up_with_time.append(array([nan]))
                self.PGD_cart_with_time.append(array([nan]))                                         
                                                                       
                                                                                     

    def get_station_magnitudes(self,A=-4.434,B=1.047,C=-0.138):
        '''
        Get PGD for every station as a function of time
        
        '''
        
        from numpy import where,c_,array,nan,log10
        
        self.station_magnitudes=[]
        
        for k in range(self.Nsta):
            
            PGD=self.PGD_cart_with_time[k]
            if len(PGD)>1:
                dist_in_km=self.hypo_distance_in_km[k]
                pgd=PGD[:,1]
                pgd=pgd*100 #to cm 'cause scaling law is in km and cm
                time=PGD[:,0]
                
                Mw=(log10(pgd)-A)/(B+C*log10(dist_in_km))
                
                #Set PGD =0 times to Mw=nan
                i=where(pgd==0)[0]
                Mw[i]=nan
                self.station_magnitudes.append(c_[time,Mw])
                                                                                                                            
            else:
                self.station_magnitudes.append(array([nan]))       
                

    def get_event_magnitudes(self,tmax=200,dt=1.0,thresh=0.02,method='mean',weighted=True):
        '''
        Get magnitudes
        
        '''
        
        from numpy import arange,zeros,where,nan,array,c_,median,average
        
        t=arange(0,tmax,dt)
        Mw=zeros(len(t))
        
        for ktime in range(len(t)):
            Mw_pgd=[]
            weights=[]
            for ksta in range(self.Nsta):
                
                pgd=self.PGD_cart_with_time[ksta]
                dist=self.hypo_distance_in_km[ksta]
                if len(pgd)>1:
                    pgd_time=pgd[:,0]
                    pgd=pgd[:,1]
                    Mw_station=self.station_magnitudes[ksta][:,1]
                    i=where(pgd_time==t[ktime])[0]
                    
                    if pgd[i]>thresh:
                        Mw_pgd.append(Mw_station[i][0])
                        weights.append(1./dist**2)
            Mw_pgd=array(Mw_pgd)
            weights=array(weights)
            
    
            
            if len(Mw_pgd)>1:
                if method=='mean':
                    if weighted==False:
                        Mw[ktime]=Mw_pgd.mean()
                    else:
                        weights=weights/max(weights)
                        Mw[ktime]=average(Mw_pgd,weights=weights)
                elif method=='median':
                    Mw[ktime]=median(Mw_pgd)
                else:
                    'ERROR: unknown method'
            else:
                Mw[ktime]=nan                                                                                                                                                                                                                           
         
        self.Mw_with_time=c_[t,Mw]                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                               


class pgd(object):
    '''
    PGD GMPE setup
    '''
    
    def __init__(self,path_to_files):
        '''
        Initalize the class
        '''
        
        from numpy import genfromtxt
        
        self.path_to_files=path_to_files
        self.catalog_file=path_to_files+'_station_info/_summary.txt'
        
        #Read event metadata
        self.event_names=genfromtxt(self.catalog_file,usecols=0,dtype='S')
        self.event_epicentral_times=genfromtxt(self.catalog_file,usecols=1,dtype='S')
        self.event_hypocenters=genfromtxt(self.catalog_file,usecols=[2,3,4])
        self.magnitudes=genfromtxt(self.catalog_file,usecols=5)
        self.Nevents=len(self.event_names)
        self.durations=genfromtxt(self.catalog_file,usecols=6)
        
        #Create event objects
        self.event=[]
        for k in range(self.Nevents):
            hypo=self.event_hypocenters[k]
            hypo_time=self.event_epicentral_times[k]
            duration=self.durations[k]
            path_to_data=self.path_to_files+self.event_names[k]
            path_to_station_file=self.path_to_files+'_station_info/'+self.event_names[k]+'.sta'
            self.event.append(Event(path_to_data,path_to_station_file,hypo,hypo_time,duration))
            
    
    
    def get_travel_times(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating phases for event '+self.event_names[k]
            self.event[k].get_travel_times()
            
            
            
    def get_SNR(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating SNR for all stations for event '+self.event_names[k]
            self.event[k].get_SNR()

    def get_PGD(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGD for all stations for event '+self.event_names[k]
            self.event[k].get_PGD() 
            
                       
    def get_PGD_with_time(self):
        '''
        Get PGDs as a function of time
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGDs with time for all stations for event '+self.event_names[k]
            self.event[k].get_PGD_with_time() 
            
    def get_station_magnitudes(self):
        '''
        Get station_mags
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGDs with time for all stations for event '+self.event_names[k]
            self.event[k].get_station_magnitudes() 
            
    def get_event_magnitudes(self,thresh=0.02,method='mean',weighted=True):
        '''
        Get event_mags
        '''
        
        for k in range(self.Nevents):
            print '... calculating Mw as a function of time for event '+self.event_names[k]
            self.event[k].get_event_magnitudes(thresh=thresh,method='mean',weighted=True)    
                      
    
    def get_hypocentral_distance(self):
        '''
        Get station to hypo distance (considers depth too)
        '''
        
        for k in range(self.Nevents):
            print '... calculating hypocentral distance for all stations for event '+self.event_names[k]
            self.event[k].get_hypocentral_distance()
            
            
    def plot_PGD(self,SNR_thresh=2.0,max_dist=1000):
        '''
        Make scatter plot
        '''
        
        from matplotlib import pyplot as plt
        import itertools
        from numpy import where,intersect1d,isnan
        
        fig=plt.figure(figsize=(12,8))
        ax=fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_xscale('log')
        
        #Colors and symbols
        marker = itertools.cycle(('v', '^', 's', 'o', 'h','D')) 
        color = itertools.cycle(('#DC143C','#228B22','#1E90FF','#4B0082','#FFD700','#32CD32','#0000CD')) 
        
        for k in range(self.Nevents):
            x=self.event[k].hypo_distance_in_km
            y=self.event[k].PGD_cart
            
            #Filter by distance
            dist=self.event[k].hypo_distance_in_km
            filter1=where(dist < max_dist)[0]
            
            #Filter by SNR
            cumulative_snr=self.event[k].SNR_north+self.event[k].SNR_east+self.event[k].SNR_up
            i=where(isnan(cumulative_snr)==1)[0]
            cumulative_snr[i]=0
            filter2=where(cumulative_snr > SNR_thresh)[0]
            
            filt=intersect1d(filter1,filter2)
            x=x[filt]
            y=y[filt]
            
            label='M'+str(self.magnitudes[k])+' ' +self.event_names[k]
            ax.scatter(x,y,marker=marker.next(),c=color.next(),s=38,label=label)
        
        ax.legend(bbox_to_anchor=(1.67, 1.0))
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('PGD (m)')
        plt.subplots_adjust(right=0.65,top=0.97,left=0.09,bottom=0.1)
        plt.show()
            
            
            


class Pwave_Event(object):
    '''
    Things related to a particular earthquakes
    '''

    def __init__(self,path_to_data,hypo,hypo_time,all_stations,all_stations_lonlat,ID):
        
        '''
        Initalize the class
        '''        
        
        from numpy import genfromtxt
        
        self.hypocenter=hypo
        self.hypo_time=hypo_time
        self.path_to_data=path_to_data
        self.all_stations=all_stations
        self.all_stations_lonlat=all_stations_lonlat
        self.ID=ID
        
        
    def get_available_stations(self):
        '''
        Which stations from the master lsit are actually in this event
        '''
        from glob import glob
        from numpy import array,where,r_
        
        self.stations=[]
        stations_list=array(glob(self.path_to_data+'*HNE*'))
        for k in range(len(stations_list)):
            sta=stations_list[k].split('/')[-1].split('.')[0]
            self.stations.append(sta)
            i=where(self.all_stations==sta)[0]
            lonlat=self.all_stations_lonlat[i,:]
            if k==0:
                self.stations_lonlat=lonlat.copy()
            else:
                self.stations_lonlat=r_[self.stations_lonlat,lonlat]
                
                
    def filter_stations_by_distance(self,dist_filter_in_km=300):
        '''
        Only keep stations within a certain distance of hypo
        '''
        
        from pyproj import Geod
        from numpy import where,ones,array
            
        proj=Geod(ellps='WGS84')
            
        #get all distances
        lat1=self.hypocenter[1]*ones(len(self.stations))
        lon1=self.hypocenter[0]*ones(len(self.stations))
        lat2=self.stations_lonlat[:,1]
        lon2=self.stations_lonlat[:,0]
        az,baz,dist_in_m=proj.inv(lon1, lat1, lon2, lat2)
        dist_in_km=dist_in_m/1000.
        
        #Filter
        i=where(dist_in_km<dist_filter_in_km)[0]
        self.filtered_stations=array(self.stations)[i]
        self.filtered_stations_lonlat=self.stations_lonlat[i,:]
        self.filtered_station_distance=dist_in_km[i]
        
    def read_waveforms(self):
        '''
        Read all waveforms
        '''
        
        from obspy import read,Stream,Trace
        
        self.waveforms_east=Stream()
        self.waveforms_north=Stream()
        self.waveforms_up=Stream()
        
        for k in range(len(self.filtered_stations)):
            east=self.path_to_data+self.filtered_stations[k]+'.HNE.sac'
            north=self.path_to_data+self.filtered_stations[k]+'.HNN.sac'
            up=self.path_to_data+self.filtered_stations[k]+'.HNZ.sac'
            
            self.waveforms_east+=read(east)
            self.waveforms_north+=read(north)
            self.waveforms_up+=read(up)
            
            
    def baseline_correct(self,fcorner=[0.075,30],N_pre_event_samples=1000,order=4,zerophase=False):
        '''
        Remove pre-event mean and bandpass filter
        '''
        
        from scipy.signal import butter,filtfilt,lfilter
        from numpy import array,mean
        
        for k in range(len(self.filtered_stations)):
            
            data1=self.waveforms_east[k].data
            data2=self.waveforms_north[k].data
            data3=self.waveforms_up[k].data
            
            #pre-event mean
            data1=data1-mean(data1[0:N_pre_event_samples])
            data2=data2-mean(data2[0:N_pre_event_samples])
            data3=data3-mean(data3[0:N_pre_event_samples])
            
            fsample=1./self.waveforms_east[k].stats.delta
            fnyquist=fsample/2
            b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
            
            if zerophase==True:
                data_filt1=filtfilt(b,a,data1)
                data_filt2=filtfilt(b,a,data2)
                data_filt3=filtfilt(b,a,data3)
            else:
                data_filt1=lfilter(b,a,data1)
                data_filt2=lfilter(b,a,data2)
                data_filt3=lfilter(b,a,data3)
            
            self.waveforms_east[k].data=data_filt1
            self.waveforms_north[k].data=data_filt2
            self.waveforms_up[k].data=data_filt3
            
            
    
    def integrate(self):
        '''
        Produce velocities and accelerations
        '''
        
        from scipy.integrate import cumtrapz
        
        self.velocities_east=self.waveforms_east.copy()
        self.displacements_east=self.waveforms_east.copy()
        self.velocities_north=self.waveforms_north.copy()
        self.displacements_north=self.waveforms_north.copy()
        self.velocities_up=self.waveforms_up.copy()
        self.displacements_up=self.waveforms_up.copy()
        
        
        for k in range(len(self.filtered_stations)):
            
            time=self.waveforms_east[k].times()
            data1=self.waveforms_east[k].data
            data2=self.waveforms_north[k].data
            data3=self.waveforms_up[k].data
            
            v1=cumtrapz(data1,time,initial=0)
            d1=cumtrapz(v1,time,initial=0)
            v2=cumtrapz(data2,time,initial=0)
            d2=cumtrapz(v2,time,initial=0)
            v3=cumtrapz(data3,time,initial=0)
            d3=cumtrapz(v3,time,initial=0)
            
            self.velocities_east[k].data=v1
            self.displacements_east[k].data=d1
            self.velocities_north[k].data=v2
            self.displacements_north[k].data=d2
            self.velocities_up[k].data=v3
            self.displacements_up[k].data=d3
            
            
    def pick(self,path_to_file=None,from_file=False,method='sta_lta',trigger_on=2.0,sta=1.0,lta=5.0):
        '''
        Pick onsets
        '''
        
        from obspy.signal.trigger import recursive_sta_lta
        from numpy import zeros,where,nan,genfromtxt,array
        
        if from_file==False:
            self.picks=zeros(len(self.filtered_stations))
        
            for k in range(len(self.filtered_stations)):
                df = self.waveforms_up[k].stats.sampling_rate
                cft = recursive_sta_lta(self.waveforms_up[k], int(sta * df), int(lta * df))
                i=where(cft>trigger_on)[0]
                if len(i)>0:
                    self.picks[k]=self.waveforms_up[k].times()[i[0]]
                else:
                    i=where(cft>trigger_on/2.)[0]
                    self.picks[k]=self.waveforms_up[k].times()[i[0]]
                    
        else:
            picks=genfromtxt(path_to_file+self.ID+'.pick',usecols=1,dtype='S')
            for k in range(len(picks)):
                if picks[k]=='nan':
                    picks[k]=float(nan)
                else:
                    picks[k]=float(picks[k]) 
                picks=array(picks).astype('float')
                self.picks=array(picks)
                
                
    def get_pgd(self,window_length=5.0):
        '''
        Get PGD from pick to window length
        '''
        
        from numpy import zeros,where,isnan,nan,sqrt
        
        self.pgd_vert=zeros(len(self.filtered_stations))
        self.pgd_horiz=zeros(len(self.filtered_stations))
        self.pgd_all=zeros(len(self.filtered_stations))
        
        for k in range(len(self.filtered_stations)):
            
            time=self.waveforms_east[k].times()
            pick=self.picks[k]
            if isnan(pick)==False:
                n=self.displacements_north[k]
                e=self.displacements_east[k]
                z=self.displacements_up[k]
                i=where((time>pick) & (time<pick+window_length))[0]
                n=n[i]
                e=e[i]
                z=z[i]
                self.pgd_vert[k]=max(abs(z))
                self.pgd_horiz[k]=max(abs(sqrt(n**2+e**2)))
                self.pgd_all[k]=max(abs(sqrt(n**2+e**2+z**2)))
            else:
                self.pgd_vert[k]=nan
                self.pgd_horiz[k]=nan
                self.pgd_all[k]=nan
            
            
    def get_pgv(self,window_length=5.0):
        '''
        Get PGD from pick to window length
        '''
        
        from numpy import zeros,where,isnan,nan,sqrt
        
        self.pgv_vert=zeros(len(self.filtered_stations))
        self.pgv_horiz=zeros(len(self.filtered_stations))
        self.pgv_all=zeros(len(self.filtered_stations))
        
        for k in range(len(self.filtered_stations)):
            
            time=self.waveforms_east[k].times()
            pick=self.picks[k]
            if isnan(pick)==False:
                n=self.displacements_north[k]
                e=self.displacements_east[k]
                z=self.displacements_up[k]
                i=where((time>pick) & (time<pick+window_length))[0]
                n=n[i]
                e=e[i]
                z=z[i]
                self.pgv_vert[k]=max(abs(z))
                self.pgv_horiz[k]=max(abs(sqrt(n**2+e**2)))
                self.pgv_all[k]=max(abs(sqrt(n**2+e**2+z**2)))
            else:
                self.pgv_vert[k]=nan
                self.pgv_horiz[k]=nan
                self.pgv_all[k]=nan
            
            
    def get_pga(self,window_length=5.0):
        '''
        Get pga from pick to window length
        '''
        
        from numpy import zeros,where,isnan,nan,sqrt
        
        self.pga_vert=zeros(len(self.filtered_stations))
        self.pga_horiz=zeros(len(self.filtered_stations))
        self.pga_all=zeros(len(self.filtered_stations))
        
        for k in range(len(self.filtered_stations)):
            
            time=self.waveforms_east[k].times()
            pick=self.picks[k]
            if isnan(pick)==False:
                n=self.displacements_north[k]
                e=self.displacements_east[k]
                z=self.displacements_up[k]
                i=where((time>pick) & (time<pick+window_length))[0]
                n=n[i]
                e=e[i]
                z=z[i]
                self.pga_vert[k]=max(abs(z))
                self.pga_horiz[k]=max(abs(sqrt(n**2+e**2)))
                self.pga_all[k]=max(abs(sqrt(n**2+e**2+z**2)))
            else:
                self.pga_vert[k]=nan
                self.pga_horiz[k]=nan
                self.pga_all[k]=nan
            
            
            
    def dump_picks(self,fout):
        '''
        Write picks to file
        '''
        
        f=open(fout,'w')
        
        for k in range(len(self.filtered_stations)):
            sta=self.filtered_stations[k]
            pick=self.picks[k]
            line='%s\t%.4f\n' % (sta,pick)
            f.write(line)
            
        f.close()
            
        



        

        
    
class pwave_pgd(object):
    '''
    PGD GMPE setup
    '''
    
    def __init__(self,path_to_files,stations,dist_filter_in_km=300):
        '''
        Initalize the class
        '''
        
        from numpy import genfromtxt
        
        print 'Loading and pre-processing data for all events'
        
        self.path_to_files=path_to_files
        self.catalog_file=path_to_files+'catalog/events.txt'
        self.stations_file=stations
        
        ##Read event metadata
        self.event_IDs=genfromtxt(self.catalog_file,usecols=0,dtype='S')
        self.event_epicentral_times=genfromtxt(self.catalog_file,usecols=1,dtype='S')
        self.event_hypocenters=genfromtxt(self.catalog_file,usecols=[2,3,4])
        self.magnitudes=genfromtxt(self.catalog_file,usecols=5)
        self.stations=genfromtxt(self.stations_file,usecols=0,dtype='S')
        self.stations_lonlat=genfromtxt(self.stations_file,usecols=[1,2])
        self.Nevents=len(self.event_IDs)
        
        #Create event objects
        self.event=[]
        for k in range(self.Nevents):
            hypo=self.event_hypocenters[k]
            hypo_time=self.event_epicentral_times[k]
            path_to_data=self.path_to_files+'sm/'+self.event_IDs[k]+'/sac/'
            ID=self.event_IDs[k]
            self.event.append(Pwave_Event(path_to_data,hypo,hypo_time,self.stations,self.stations_lonlat,ID))
            self.event[k].get_available_stations()
            self.event[k].filter_stations_by_distance(dist_filter_in_km=dist_filter_in_km)
            print '... reading waveforms for '+self.event_IDs[k]
            self.event[k].read_waveforms()
            print '... baseline correction for '+self.event_IDs[k]
            self.event[k].baseline_correct()
            print '... integrating waveforms for '+self.event_IDs[k]
            self.event[k].integrate()
            print '... picking for '+self.event_IDs[k]
            pick_path=u'/Users/dmelgar/pwave_scaling/picks/clean/'
            self.event[k].pick(path_to_file=pick_path,from_file=True)
            #self.event[k].pick()
            
            
    
    
    def get_travel_times(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating phases for event '+self.event_names[k]
            self.event[k].get_travel_times()
            
            
            
    def get_SNR(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating SNR for all stations for event '+self.event_names[k]
            self.event[k].get_SNR()

    def get_PGD(self):
        '''
        Get travel times for all events
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGD for all stations for event '+self.event_names[k]
            self.event[k].get_PGD() 
            
                       
    def get_PGD_with_time(self):
        '''
        Get PGDs as a function of time
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGDs with time for all stations for event '+self.event_names[k]
            self.event[k].get_PGD_with_time() 
            
    def get_station_magnitudes(self):
        '''
        Get station_mags
        '''
        
        for k in range(self.Nevents):
            print '... calculating PGDs with time for all stations for event '+self.event_names[k]
            self.event[k].get_station_magnitudes() 
            
    def get_event_magnitudes(self,thresh=0.02,method='mean',weighted=True):
        '''
        Get event_mags
        '''
        
        for k in range(self.Nevents):
            print '... calculating Mw as a function of time for event '+self.event_names[k]
            self.event[k].get_event_magnitudes(thresh=thresh,method='mean',weighted=True)    
                      
    
    def get_hypocentral_distance(self):
        '''
        Get station to hypo distance (considers depth too)
        '''
        
        for k in range(self.Nevents):
            print '... calculating hypocentral distance for all stations for event '+self.event_names[k]
            self.event[k].get_hypocentral_distance()
            
            
    def plot_PGD(self,SNR_thresh=2.0,max_dist=1000):
        '''
        Make scatter plot
        '''
        
        from matplotlib import pyplot as plt
        import itertools
        from numpy import where,intersect1d,isnan
        
        fig=plt.figure(figsize=(12,8))
        ax=fig.add_subplot(111)
        ax.set_yscale('log')
        ax.set_xscale('log')
        
        #Colors and symbols
        marker = itertools.cycle(('v', '^', 's', 'o', 'h','D')) 
        color = itertools.cycle(('#DC143C','#228B22','#1E90FF','#4B0082','#FFD700','#32CD32','#0000CD')) 
        
        for k in range(self.Nevents):
            x=self.event[k].hypo_distance_in_km
            y=self.event[k].PGD_cart
            
            #Filter by distance
            dist=self.event[k].hypo_distance_in_km
            filter1=where(dist < max_dist)[0]
            
            #Filter by SNR
            cumulative_snr=self.event[k].SNR_north+self.event[k].SNR_east+self.event[k].SNR_up
            i=where(isnan(cumulative_snr)==1)[0]
            cumulative_snr[i]=0
            filter2=where(cumulative_snr > SNR_thresh)[0]
            
            filt=intersect1d(filter1,filter2)
            x=x[filt]
            y=y[filt]
            
            label='M'+str(self.magnitudes[k])+' ' +self.event_names[k]
            ax.scatter(x,y,marker=marker.next(),c=color.next(),s=38,label=label)
        
        ax.legend(bbox_to_anchor=(1.67, 1.0))
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('PGD (m)')
        plt.subplots_adjust(right=0.65,top=0.97,left=0.09,bottom=0.1)
        plt.show()    