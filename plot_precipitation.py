# Plot Precipitation (Python version)
import os
import numpy as np
#import math
import matplotlib
import glob
#matplotlib.use('Agg')
import matplotlib.pyplot as plt


#---------- Find which files we have to work with: ------------

dataDir = 'Out_45deg';

phi_files = glob.glob(dataDir + '/phi_*')
f = [f.split('/phi_')[-1] for f in phi_files]
L_NS_tuples = [l.split('_') for l in f]  # [L shell, N or S]

L_NS_tuples.sort(key=lambda val: val[0])
L_vec = np.unique(np.array([float(l[0]) for l in L_NS_tuples]))
#print L_vec

# Model constants:

ev2joule = 1.60217657*pow(10,-19) # Joules / ev
joule2millierg = pow(10,10)
# Constants for energy scaling:
# (Copied from readOne file)
DE_EXP   =  0.003
E_EXP_BOT=  1.477
E_EXP_TOP=  7.477  
RES_DT   =  0.02  
RES_FINT =  5.0
EA_SPLIT =  1
MULT     =  2.0
NUM_E    =  2200

# Energies corresponding to energy bins (columns) of data in phi file
# Energies step exponentially and are in [eV].  
E_EXP = E_EXP_BOT + np.arange(1,NUM_E+1)*DE_EXP
E = pow(10,E_EXP)    #Energy Values (eV)
DE = np.gradient(E)   #Energy Step Sizes

# time axis:
t = np.arange(RES_DT,RES_FINT,RES_DT)

N = np.zeros((t.size, NUM_E, L_vec.size))
S = np.zeros((t.size, NUM_E, L_vec.size))

for idx, L in enumerate(L_vec):
    #print idx, L
    for NS in ['N','S']:
        fn = '%s/phi_%g_%s' % (dataDir, L, NS)
     
        print(fn)

        d = np.fromfile(fn,np.float32).reshape(t.size,NUM_E,order='F')
        if NS=='N':
            N[:,:,idx] = d
        elif NS=='S':
            S[:,:,idx] = d

#------------- Huzzah! Data loaded ----------------------
# N and S are now 3d matrices (time x energy x L-shell)
# with units of [counts / (cm^2 keV s)]
# -------------------------------------------------------

E_scaled = ev2joule*joule2millierg*E
N_scaled = np.zeros((t.size, NUM_E, L_vec.size))
S_scaled = np.zeros((t.size, NUM_E, L_vec.size))

E_repped = np.array(np.tile(E_scaled, (np.shape(N)[0],1)))

for idx in range(np.shape(N)[2]):
    N_scaled[:,:,idx] = N[:,:,idx]*E_repped
    S_scaled[:,:,idx] = S[:,:,idx]*E_repped

N_sum = np.sum(N_scaled,axis=1) + pow(10,-10)
S_sum = np.sum(S_scaled,axis=1) + pow(10,-10)
# fig = plt.figure()
# plt.imshow(np.log10(E_repped))
# plt.show()

N_log = np.log10(N_scaled) # + np.log10(ev2joule*joule2millierg)
S_log = np.log10(S_scaled) # + np.log10(ev2joule*joule2millierg)
N_log[N==0] = -100
S_log[S==0] = -100
# ------------ Total flux at each L-shell (Tile plot, time vs energy)
H = int(np.floor(np.sqrt(L_vec.size)))
W = int(np.ceil(np.sqrt(L_vec.size)))
print("H, W: ",H,W)
fig, axvec = plt.subplots(W,H, sharex=True, sharey=True)
print(np.shape(axvec))
axvec = axvec.reshape(W*H,1)
#axvec = axvec(:)
for idx, L in enumerate(L_vec):
    print("idx is ", idx)
    sp = axvec[idx][0] # Current subplot
    im = sp.imshow(N_log[:,:,idx].T,aspect='auto',interpolation='none')
    plt.tight_layout()
    im.set_clim(-4, 0)
    sp.invert_yaxis()
    sp.set_title(L)
    print(sp.get_xticklabels())
    #plt.yticks(E_scaled)
    plt.setp(sp.get_xticklabels(), visible=False)
    plt.setp(sp.get_yticklabels(), visible=False)

#    xlim([0 5]);
#    caxis([-10 0]);
#    title(L_vec(l));
#    set(axvec{l},'YTick',[]);
# end
plt.suptitle('Flux per L-shell');
axvec = axvec.reshape(W,H)

plt.setp([a.get_xticklabels() for a in axvec[H,:]], visible=True)
plt.setp([a.get_yticklabels() for a in axvec[:,0]], visible=True)

fig2, av2 = plt.subplots(2,1, sharex=True, sharey=True)
for idx in (0,1):
    sp = av2[idx] # Current subplot
    im = sp.imshow(np.log10(N_sum).T,aspect='auto',interpolation='none')
    sp.invert_yaxis()
    sp.set_title(['N','S'][idx])
plt.show()


