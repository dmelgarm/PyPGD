from PyPGD import pgd
import cPickle as pickle

path_to_files=u'/Users/dmelgar/PGD_gmpe/'
file_out=u'/Users/dmelgar/PGD_gmpe/all_events_pgd.pkl'


with open(file_out, 'rb') as input_file:
    p = pickle.load(input_file)
        
