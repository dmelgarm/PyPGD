    
class Event(object):
    '''
    Things related to a particular earthquakes
    '''

    def __init__(self,path_to_data,path_to_station_file,hypocenter,hypo_time):
        
        '''
        Initalize the class
        '''        
        
        from numpy import genfromtxt
        
        self.hypocenter=hypocenter
        self.hypo_time=hypo_time
        self.path_to_data=path_to_data
        self.path_to_station_file=path_to_station_file
        self.station_locations=genfromtxt(self.path_to_station_file,usecols=[1,2])
        self.station_names=genfromtxt(self.path_to_station_file,usecols=0,dtype='S')
        
    def get_travel_times(self):
        '''
        Get direct P and S to each station
        '''
        
        from obspy.taup import TauPyModel
        from numpy import float64,ones
        from obspy.geodetics import locations2degrees
        from datetime import timedelta
        
        Nsta=len(self.station_names)

        self.ptime=1e6*ones(Nsta)
        self.stime=1e6*ones(Nsta)        
                        
        for k in range(Nsta):
        
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
        
        #Create event objects
        self.event=[]
        for k in range(self.Nevents):
            hypo=self.event_hypocenters[k]
            hypo_time=self.event_epicentral_times
            path_to_data=self.path_to_files+self.event_names[k]
            path_to_station_file=self.path_to_files+'_station_info/'+self.event_names[k]+'.sta'
            self.event.append(Event(path_to_data,path_to_station_file,hypo,hypo_time))

        

        
        