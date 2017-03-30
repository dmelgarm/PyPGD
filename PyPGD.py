    
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
        
        self.SNR_north=zeros(self.Nsta)  
        self.SNR_east=zeros(self.Nsta)   
        self.SNR_up=zeros(self.Nsta) 
        
        hypo_time=UTCDateTime(self.hypo_time)
        
        for k in range(self.Nsta):
            path=self.path_to_data
            
            ptime=self.ptime[k]

            tp=hypo_time+timedelta(seconds=ptime)  
            t0=tp-timedelta(seconds=t_prior)
            tfinal=tp+timedelta(seconds=source_durations*self.duration)          
                                    
            #read waveforms
            n=read(path+self.station_names[k]+'.LXN.sac')
            e=read(path+self.station_names[k]+'.LXE.sac')
            z=read(path+self.station_names[k]+'.LXZ.sac')
            
            #Get pre-event portions (from t0 to tp)
            n_pre=n.copy()
            e_pre=e.copy()
            z_pre=z.copy()
            
            n_pre[0].trim(starttime=t0,endtime=tp)
            e_pre[0].trim(starttime=t0,endtime=tp)
            z_pre[0].trim(starttime=t0,endtime=tp)
            
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
            
            n[0].data=n[0].data-mean_n
            e[0].data=e[0].data-mean_e
            z[0].data=z[0].data-mean_z
            
            #get power
            Pn=sum(n[0].data**2)/n[0].stats.npts
            Pe=sum(e[0].data**2)/e[0].stats.npts
            Pz=sum(z[0].data**2)/z[0].stats.npts
            
            #And finaly SNR
            self.SNR_north[k]=Pn/Pn_pre
            self.SNR_east[k]=Pe/Pe_pre
            self.SNR_up[k]=Pz/Pz_pre
               


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

        

        
        