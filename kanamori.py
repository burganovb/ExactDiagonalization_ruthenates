import numpy as np

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
    print statediff, np.sum(np.abs(statediff)), np.sum(statediff)
    if np.sum(np.abs(statediff))!=4 or np.sum(statediff)!=0:
        return 0.
    n_sites = state1.shape[0]/6
    diffind = np.argwhere(statediff)
    print n_sites, diffind
    print diffind[0],diffind[2]-3, diffind[1],diffind[3]-3
    print diffind[0]/6, diffind[-1]/6
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

#print Kanamori_offdiag_energy(np.array([0, 1, 0, 1, 0, 1]), np.array([1, 0, 0, 0, 1, 1]), 1.)
#print Kanamori_offdiag_energy(np.array([0, 1, 1, 1, 0, 1]), np.array([1, 0, 1, 0, 1, 1]), 1.)
print Kanamori_offdiag_energy(np.array([0,1, 1,0, 1,0, 1,0, 0,0, 1,0]), np.array([1,0, 0,0, 1,0, 0,0, 1,0, 1,0]), 1.)
#print Kanamori_offdiag_energy(np.array([1, 1, 0, 0, 0, 1]), np.array([0, 0, 0, 1, 1, 1]), 1.)
