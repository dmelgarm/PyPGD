from PyPGD import pgd

path_to_files=u'/Users/dmelgar/PGD_gmpe/'
p=pgd(path_to_files)


# Get hypocentral distance too DOESN'T EXIST YET
p.get_hypocentral_distance()

#get travel times THIS IS NOT WORKING FOR ALL EVENTS tohoku lon lat columns are transposed
p.get_travel_times()

#get SNR
p.get_SNR()