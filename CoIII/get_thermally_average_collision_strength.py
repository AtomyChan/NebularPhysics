import os
import numpy as np
from astropy import constants
from astropy import units

def get_collision_strength(i,j, datadir='/Users/chenping/Ping/NebularPhysics'):
	

	fname = 'Storey2016/onlinedata/data/OMEGA_%s_%s_CoIII.dat'%(i,j)
	datafile = os.path.join(datadir, fname)
	
	data = np.loadtxt(datafile)

	return data



def thermal_average(i,j,T):

	csdata = get_collision_strength(i,j)
	kT_eV = (constants.k_B*T*units.K).to('eV').value

	TavCS = 0

	N = len(csdata)

	for i,ci in enumerate(csdata):
		if i<(N-1):	
	
			csi = ci[1]
			Ei = ci[0]
			cip1 = csdata[i+1]
			deltaE = cip1[0] - ci[0]
	
			inti = csi*np.exp(-Ei/kT_eV)*deltaE/kT_eV
	
			TavCS = TavCS + inti

	return TavCS


if __name__ == '__main__':
	import sys
	
	i = sys.argv[1]
	j = sys.argv[2]
	T = float(sys.argv[3])

	TavCS = thermal_average(i,j, T)
	print TavCS
