import numpy as np
from scipy.special import factorial
from numpy.linalg import eigh
import matplotlib.pyplot as plt
import blockdiag

def generate_orb_states(n_el, n_sites):
    if n_el == 0:
        return np.zeros((1, n_sites))
    elif n_el == n_sites:
        return np.ones((1, n_sites))
    else:
        states = np.zeros((0, n_sites))
        for last_pos in range(n_el-1, n_sites):
            tmp = generate_orb_states(n_el-1, last_pos)
            ones_arr = np.ones((tmp.shape[0], min(1,n_sites-last_pos)))
            zeros_arr = np.zeros((tmp.shape[0], max(0,n_sites-1-last_pos)))
            tmp1 = np.concatenate((np.concatenate((tmp,ones_arr),axis=1),zeros_arr),axis=1)
            states = np.concatenate((states, tmp1), axis=0)
        return states

def generate_all_states(n_orb_spin, n_sites):
    # output is #rows=#states
    # first n_sites columns = occupancy of xy-up electrons on diff-t sites
    # next n_sites columns = occupancy of xz-up electrons on diff-t sites, etc.
    xyup = generate_orb_states(n_orb_spin[0], n_sites)
    xzup = generate_orb_states(n_orb_spin[1], n_sites)
    yzup = generate_orb_states(n_orb_spin[2], n_sites)
    xydn = generate_orb_states(n_orb_spin[3], n_sites)
    xzdn = generate_orb_states(n_orb_spin[4], n_sites)
    yzdn = generate_orb_states(n_orb_spin[5], n_sites)
    n_states = xyup.shape[0]*xzup.shape[0]*yzup.shape[0]*xydn.shape[0]*xzdn.shape[0]*yzdn.shape[0]
    states = np.zeros((n_states, 6*n_sites))
    k = 0
    for k0 in range(xyup.shape[0]):
        for k1 in range(xzup.shape[0]):
            for k2 in range(yzup.shape[0]):
                for k3 in range(xydn.shape[0]):
                    for k4 in range(xzdn.shape[0]):
                        for k5 in range(yzdn.shape[0]):
                            states[k, 0:n_sites] = xyup[k0]
                            states[k, n_sites:2*n_sites] = xzup[k1]
                            states[k, 2*n_sites:3*n_sites] = yzup[k2]
                            states[k, 3*n_sites:4*n_sites] = xydn[k3]
                            states[k, 4*n_sites:5*n_sites] = xzdn[k4]
                            states[k, 5*n_sites:6*n_sites] = yzdn[k5]
                            k+=1
    return states

def onsite_energy(state, U, J):
    n_sites = state.shape[0]/6
    #print n_sites
    energy = 0
    for k in range(n_sites):
        n_el = np.sum(state[k::n_sites])
        n_up = np.sum(state[k:3*n_sites:n_sites])
        n_dn = np.sum(state[3*n_sites+k:6*n_sites:n_sites])
        #print n_el, 0.5*U*n_el*(n_el-1)
        #energy += 0.5*U*n_el*(n_el-1) + J*(n_up*n_dn - 0.5*n_up*(n_up-1) - 0.5*n_dn*(n_dn-1))
        energy += (U-1.5*J)*0.5*n_el*(n_el-1) - 0.25*J*(n_up-n_dn)**2
    return energy
    
# PRB 80, 235114
def Kanamori_diag_energy(state, U, J):
    n_sites = state.shape[0]/6
    energy = 0
    for k in range(n_sites):
        siteoccup = state[k::n_sites]
        energy += U*(siteoccup[0]*siteoccup[3]+siteoccup[1]*siteoccup[4]+siteoccup[2]*siteoccup[5]) #U
        energy += (U-3*J)*(siteoccup[0]*(siteoccup[1]+siteoccup[2])+siteoccup[1]*siteoccup[2]) #up's
        energy += (U-3*J)*(siteoccup[3]*(siteoccup[4]+siteoccup[5])+siteoccup[4]*siteoccup[5]) #dn's
        energy += (U-2*J)*(siteoccup[0]*(siteoccup[4]+siteoccup[5])+siteoccup[1]*(siteoccup[3]+siteoccup[5])+siteoccup[2]*(siteoccup[3]+siteoccup[4])) #up-dn
    return energy

def Kanamori_offdiag_energy(state1, state2, J):
    statediff = state1 - state2
    if np.sum(np.abs(statediff))!=4 or np.sum(statediff)!=0:
        return 0.
    n_sites = state1.shape[0]/6
    diffind = np.argwhere(statediff)
    siteind = np.mod(diffind[0],6)
    orb_occup_diff = statediff[siteind::n_sites]
    
    if np.sum(orb_occup_diff[0:3])!=0 or np.sum(np.abs(orb_occup_diff[0:3]))!=2 or np.sum(orb_occup_diff[3:6])!=0 or np.sum(np.abs(orb_occup_diff[3:6]))!=2:
        return 0. #
    if all(orb_occup_diff[0:3]==-orb_occup_diff[3:6]):
        return -J # spin exchange
    elif all(orb_occup_diff[0:3]==orb_occup_diff[3:6]):
        return -J #pairing
    else:
        return 0.
        
def nn_hopping_term(state1, state2, site_coords): #nearest neighbor hopping
    statediff = state1 - state2
    if np.sum(np.abs(statediff))!=2:
        return 0.
    else:
        diff1 = np.argwhere(statediff<0)
        diff2 = np.argwhere(statediff>0)
        if diff1/site_coords.shape[0] != diff2/site_coords.shape[0]: #changes in different orbitals (this should never occur at this point)
            print 'smth wrong with states'
            return 0.
        orbital = diff1/site_coords.shape[0] # from 0 to 5: xy-up,xz-up,yz-up,xy-dn,xz-dn,yz-dn
        coord1 = np.mod(diff1,site_coords.shape[0])
        coord2 = np.mod(diff2,site_coords.shape[0])
        #print orbital, coord1, coord2
        if (np.abs(site_coords[coord1] - site_coords[coord2])==np.array([1,0,0])).all() and orbital in [0,1,3,4]:
            return -1.
        elif (np.abs(site_coords[coord1] - site_coords[coord2])==np.array([0,1,0])).all() and orbital in [0,2,3,5]:
            return -1.
        elif (np.abs(site_coords[coord1] - site_coords[coord2])==np.array([0,0,1])).all() and orbital in [1,2,4,5]:
            return -1.
        else:
            return 0.

def calc_eigenstates(site_coords, n_orb_spin, U, J, tvals):
    if site_coords.shape[1]==2:
        site_coords = np.concatenate((site_coords, np.zeros((site_coords.shape[0],1))),axis=1)
    n_sites = site_coords.shape[0]
    states = generate_all_states(n_orb_spin, n_sites)
    n_states = states.shape[0]
    
    print n_states
    #print states
    
#    Hamilt_onsite = np.zeros((n_states, n_states))
    Hamilt_kanamori = np.zeros((n_states, n_states))
    Hamilt_hop = np.zeros((n_states, n_states))
#    Hamilt_offdiag = np.zeros((n_states, n_states))
    for k in range(n_states):
        #Hamilt_onsite[k,k] = onsite_energy(states[k,:], U, J)
        Hamilt_kanamori[k,k] = Kanamori_diag_energy(states[k,:], U, J)
    
    #print np.diag(Hamilt_kanamori)#-min(np.diag(Hamilt_kanamori))
    
    for k0 in range(n_states):
        for k1 in range(k0):
            hop = nn_hopping_term(states[k0,:], states[k1,:], site_coords)
            kanamoffdiag = Kanamori_offdiag_energy(states[k0,:], states[k1,:], J)
            Hamilt_kanamori[k0,k1] = kanamoffdiag
            Hamilt_kanamori[k1,k0] = kanamoffdiag
            Hamilt_hop[k0,k1] = hop
            Hamilt_hop[k1,k0] = hop
    eigvals = np.zeros((tvals.size, n_states))

    #print Hamilt_onsite
    #print Hamilt_kanamori
    #print Hamilt_hop
    print Hamilt_kanamori.shape
    k=0
    
    #####
    Hamilt = Hamilt_kanamori + 0.1*Hamilt_hop
    return Hamilt
    #####
    
    for t in tvals:
        Hamilt = Hamilt_kanamori + t*Hamilt_hop
        [w, v] = eigh(Hamilt)
        eigvals[k] = w
        k+=1
    
    return eigvals

#========== main start ===================

n_sites = 2
site_coords = np.array([[0,0], [1,0]])
n_orb_spin1 = np.array([2, 2, 2, 1, 1, 0]) # up-xy,xz,yz, dn-xy,xz,yz
n_orb_spin2 = np.array([1, 1, 2, 1, 1, 2]) # up-xy,xz,yz, dn-xy,xz,yz
#n_orb_spin2 = np.array([2, 1, 1, 1, 2, 1]) # up-xy,xz,yz, dn-xy,xz,yz
#n_orb_spin2 = np.array([2, 1, 1, 1, 1, 2]) # up-xy,xz,yz, dn-xy,xz,yz

'''
n_sites = 4
site_coords = np.array([[0,0], [1,0], [0,1], [1,1]])
n_orb_spin1 = np.array([4, 4, 4, 1, 2, 1]) # up-xy,xz,yz, dn-xy,xz,yz
n_orb_spin2 = np.array([2, 3, 3, 2, 3, 3]) # up-xy,xz,yz, dn-xy,xz,yz
'''
'''
n_sites = 3
site_coords = np.array([[0,0], [1,0], [2,0]])
n_orb_spin1 = np.array([3, 3, 3, 1, 1, 1]) # 27
#n_orb_spin2 = np.array([2, 2, 1, 2, 2, 3]) # 243
n_orb_spin2 = np.array([2, 2, 1, 2, 3, 2]) # 243
'''
n_sites = 4
site_coords = np.array([[0,0], [1,0], [2,0], [3,0]])
n_orb_spin1 = np.array([4, 4, 4, 1, 2, 1]) # 96
n_orb_spin2 = np.array([2, 2, 4, 2, 2, 4]) # 9216
#n_orb_spin1 = np.array([2, 3, 3, 2, 3, 3]) # 
#n_orb_spin2 = np.array([3, 3, 2, 3, 3, 2]) # 


#tvals = np.arange(0.0,1.6,0.05)
tvals = np.arange(0.0,0.6,0.1)
#eigvals1 = calc_eigenstates(site_coords, n_orb_spin1, 2.0, 0., tvals)
#eigvals2 = calc_eigenstates(site_coords, n_orb_spin2, 2.0, 0., tvals)

eigvals1 = calc_eigenstates(site_coords, n_orb_spin1, 2.0, 0.4, tvals)
#print eigvals1.shape
#blockdiag.blockdiag2blocks(eigvals1)

eigvals2 = calc_eigenstates(site_coords, n_orb_spin2, 2.0, 0.4, tvals)


#print eigvals1.shape
#print eigvals2[:,0:eigvals1.shape[1]].shape
'''
plt.figure()

eigvals = eigvals1-eigvals2[:,0:eigvals1.shape[1]]
plt.plot(tvals, eigvals1[:,0]-eigvals2[:,0])
#plt.plot(tvals, eigvals1[:,0], 'r')
#plt.plot(tvals, eigvals2[:,0], 'b')

plt.figure()
for k in range(eigvals1.shape[1]):
    plt.plot(tvals**2, eigvals1[:,k],'r')
for k in range(eigvals2.shape[1]):
    plt.plot(tvals**2, eigvals2[:,k],'b')
plt.figure()
for k in range(4):
    plt.plot(tvals**2, eigvals1[:,k],'r')
    plt.plot(tvals**2, eigvals2[:,k],'b')

plt.show()
'''
#print eigvals1[:,0]
#print eigvals2[:,0]