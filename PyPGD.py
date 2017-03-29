class pgd:
    '''
    Read an manipulate the USGS NEIC finite fault database
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
        
        
        