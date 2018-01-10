from PyPGD import pwave_pgd
import cPickle as pickle
from matplotlib import pyplot as plt
from numpy import reshape
from matplotlib.ticker import MultipleLocator

path_to_files=u'/Users/dmelgar/pwave_scaling/'
stations=u'/Users/dmelgar/pwave_scaling/catalog/knet.sta'
file_out=u'/Users/dmelgar/pwave_scaling/all_events_pgd.pkl'

#Do everything or load from fie?
load_from_file=False

if load_from_file==False:

    p=pwave_pgd(path_to_files,stations)
    
    #Save the current object
    with open(file_out, 'wb') as output:
        pickle.dump(p, output, pickle.HIGHEST_PROTOCOL)
        
        
    #dump picks
    for k in range(len(p.event_IDs)):
        fout='/Users/dmelgar/pwave_scaling/picks/auto/'+p.event_IDs[k]+'.pick'
        p.event[k].dump_picks(fout)
        
else: #load from file
    with open(file_out, 'rb') as input_file:
        p = pickle.load(input_file)
        



#make scatter plot
e1=p.event[2]
e2=p.event[5]
#e3=p.event[5]

lab1='M9.0'
lab2='M8.3'

#windows=[5,10,15,20,25,30,35,40,45]
windows=[10,20,30,45,60,75,90,105,120]


xmajorLocator = MultipleLocator(xtick)


s=8
xl=[50,310]
yl=[0.01,45]
fig, axarr = plt.subplots(3, 3,figsize=(16,16))
axarr=reshape(axarr,9,1)

for k in range(len(axarr)):
    
    ax=axarr[k]
    ax.set_yscale('log')
    #ax.set_xscale('log')
    e1.get_pgd(window_length=windows[k])
    e2.get_pgd(window_length=windows[k])
    #e3.get_pgd(window_length=windows[k])
    ax.scatter(e1.filtered_station_distance,e1.pgd_vert*100,s=s,label=lab1)
    ax.scatter(e2.filtered_station_distance,e2.pgd_vert*100,s=s,label=lab2)
    #ax.xaxis.set_ticklabels(['60','100','200','300'])
    
    #ax.scatter(e3.filtered_station_distance,e3.pgd*100,s=s,label='M8.3')
    ax.legend(loc=3)
    ax.set_title(str(windows[k])+'s',fontsize=14)
    if k==2 or k==5 or k==8:
        ax.set_xlabel('d (km)')
    if k<3:
        ax.set_ylabel('pgd (cm)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    
plt.subplots_adjust(hspace=0.3,left=0.1,right=0.98,bottom=0.1,top=0.95)
    
    
    
    
s=8
xl=[50,310]
yl=[0.01,45]
fig, axarr = plt.subplots(3, 3,figsize=(16,16))
axarr=reshape(axarr,9,1)

for k in range(len(axarr)):
    
    ax=axarr[k]
    ax.set_yscale('log')
    ax.set_xscale('log')
    e1.get_pgd(window_length=windows[k])
    e2.get_pgd(window_length=windows[k])
    #e3.get_pgv(window_length=windows[k])
    ax.scatter(e1.filtered_station_distance,e1.pgd_horiz*100,s=s,label=lab1)
    ax.scatter(e2.filtered_station_distance,e2.pgd_horiz*100,s=s,label=lab2)
    #ax.scatter(e3.filtered_station_distance,e3.pgv*100,s=s,label='M7.3')
    ax.legend(loc=3)
    ax.set_title('pgv at '+str(windows[k])+'s')
    ax.set_xlabel('d (km)')
    ax.set_ylabel('pgv (cm/s)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    
    
    
    
    
s=8
xl=[50,310]
yl=[0.01,45]
fig, axarr = plt.subplots(3, 3,figsize=(16,16))
axarr=reshape(axarr,9,1)

for k in range(len(axarr)):
    
    ax=axarr[k]
    ax.set_yscale('log')
    ax.set_xscale('log')
    e1.get_pgd(window_length=windows[k])
    e2.get_pgd(window_length=windows[k])
    #e3.get_pga(window_length=windows[k])
    ax.scatter(e1.filtered_station_distance,e1.pgd_all*100,s=s,label=lab1)
    ax.scatter(e2.filtered_station_distance,e2.pgd_all*100,s=s,label=lab2)
    #ax.scatter(e3.filtered_station_distance,e3.pga*100,s=s,label='M7.3')
    ax.legend(loc=3)
    ax.set_title('pga at '+str(windows[k])+'s')
    ax.set_xlabel('d (km)')
    ax.set_ylabel('pga (cm/s/s)')
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    
plt.show()