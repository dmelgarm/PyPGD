from PyPGD import pgd
import cPickle as pickle

path_to_files=u'/Users/dmelgar/PGD_gmpe/'
file_out=u'/Users/dmelgar/PGD_gmpe/all_events_pgd.pkl'

#Do everything or load from fie?
load_from_file=True

if load_from_file==False:

    p=pgd(path_to_files)
    
    # Get hypocentral distance too DOESN'T EXIST YET
    p.get_hypocentral_distance()
    
    #get travel times THIS IS NOT WORKING FOR ALL EVENTS tohoku lon lat columns are transposed
    p.get_travel_times()
    
    #get SNR
    p.get_SNR()
    
    #get PGD
    p.get_PGD()
    
    #get PGD with time
    p.get_PGD_with_time()
    
    #get PGD with time
    p.get_station_magnitudes()
    
    #get PGD with time
    p.get_event_magnitudes(thresh=0.025,method='mean',weighted=True)
    
    #plot
    p.plot_PGD(SNR_thresh=3.0)
    
    #Save the current object
    with open(file_out, 'wb') as output:
        pickle.dump(p, output, pickle.HIGHEST_PROTOCOL)
        
else: #load from file
    with open(file_out, 'rb') as input_file:
        p = pickle.load(input_file)
        
    #plot
    p.plot_PGD(SNR_thresh=3.0)