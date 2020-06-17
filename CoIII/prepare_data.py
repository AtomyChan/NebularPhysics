import numpy as np
from astropy.table import Table
from astropy import units
from astropy import constants

def reformat_table_transition_probablities(infile):
	
	jlist = []
	ilist = []
	Alist = []	

	alllines = open(infile).readlines()
	Nlines = len(alllines)	

	for k,line in enumerate(alllines):
		linesegs = line.split('&')
		j = int(linesegs[0].replace(' ',''))
		i = int(linesegs[1].replace(' ',''))
		A_str = linesegs[2].replace(' ','')
		A = float(A_str.split('(')[0])*10**(float(A_str[-3:-1]))
		jlist.append(j)
		ilist.append(i)
		Alist.append(A)
		if k<(Nlines-1):		
			j = int(linesegs[7].replace(' ',''))
			i = int(linesegs[8].replace(' ',''))
			A_str = linesegs[9].replace(' ','')
			A = float(A_str.split('(')[0])*10**(float(A_str[-3:-1]))
			
			jlist.append(j)
			ilist.append(i)
			Alist.append(A)

	data = [jlist, ilist, Alist]
	colnames = ['j','i','Aji']
	datatable = Table(data, names = colnames)
	
	#print datatable

	return datatable

def get_transition_probability(j,i, Atablefile='/Users/chenping/Ping/NebularPhysics/CoIII/Atable_Storey2016.txt'):
	'''
	The Avalue table is from Storey et al. 2016 table  6
	
	INPUTS:
		j: upper level
		i: lower lvel
		
	'''
	if j<=i:
		raise ValueError("This is spontaneous transition from upper level to lower level. j>i is required")
	
	Atable = Table.read(Atablefile, format='ascii.fixed_width')
	ji = (Atable['j']==j)*(Atable['i']==i)
	
	if np.sum(ji):
		A = Atable[ji]['Aji'].data[0]
	else:
		A = 0 #if not found in the table then return 0; (only transition probabilities that are at least 1 percent of the total probability from a given upper level are listed in the table)
	return A


def get_Tave_collision_strength(i,j,T,CStablefile='/Users/chenping/Ping/NebularPhysics/CoIII/Storey2016/arxiv/Tave_collision_strength.txt'):
	'''
	Data from Storey et al. 2017 table 7
	
	'''
	log10T = np.array([2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4])
	Tline = 10**log10T		
	data = np.loadtxt(CStablefile)
	
	if j <= i:
		raise ValueError("This is collision strength of excitatoion. j>i is required")

	CSline = data[(data[:,0]==i)*(data[:,1]==j)][0][2:]
	#print CSline
	TCSinterp = np.interp(T, Tline, CSline)
	
	return TCSinterp
	

def get_transition_energy(i,j, ELfile='/Users/chenping/Ping/NebularPhysics/CoIII/Storey2016/arxiv/energy_levels_lowest15.txt', eunit='wn'):
	'''
	Data from Storey2016/arxiv/energy_levels_lowest15.txt

	The energy in the table is in cm^(-1)
	
	INPUTS:
		j:
		i:
		ELfile:
		outunit: 'WN' for wave number in cm^{-1}; 'eV' in electron Volt
	'''
	ELs = Table.read(ELfile, format='ascii.commented_header')
	iE = ELs[ELs['index']==i]['Eexp'].data[0]
	jE = ELs[ELs['index']==j]['Eexp'].data[0]

	deltaE = jE - iE

	if eunit == 'eV':
		deltaE = (constants.h*constants.c*deltaE*(1/units.cm)).to('eV').value

	return deltaE


def get_statistical_weight(i, ELfile='/Users/chenping/Ping/NebularPhysics/CoIII/Storey2016/arxiv/energy_levels_lowest15.txt'):
	'''
	statistical weight wi = 2*J+1 where is J is total angular momentum 
	'''
	ELs = Table.read(ELfile, format='ascii.commented_header')

	Jn = ELs[ELs['index']==i]['Jn'].data[0]
	Jd = ELs[ELs['index']==i]['Jd'].data[0]
	
	wi = float(Jn)/float(Jd)*2+1
	
	return wi



if __name__ == "__main__":

	import os

	cdir = os.path.dirname(os.path.realpath(__file__))
	transition_probility_file = os.path.join(cdir, 'Storey2016/arxiv/transition_probability.tex')

	#Atable = reformat_table_transition_probablities(transition_probility_file)
	Atablefile = 'Atable_Storey2016.txt'
	#Atable.write(Atablefile, format='ascii.fixed_width')

	j  = 8
	i  = 1
	Avalue = get_transition_probability(j,i, Atablefile=Atablefile)
	print Avalue



 	Tave_CS_file = os.path.join(cdir, 'Storey2016/arxiv/Tave_collision_strength.txt')
	
	T = 12000
	i = 1
	j = 8
	TaveCS = get_Tave_collision_strength(i,j,T,CStablefile=Tave_CS_file)
	print TaveCS


	ELfile = os.path.join(cdir, 'Storey2016/arxiv/energy_levels_lowest15.txt')
	j = 9
	i = 2
	deltaE = get_transition_energy(i,j, ELfile=ELfile, eunit='eV')
	print deltaE
	
	w8 = get_statistical_weight(8)
	print w8
