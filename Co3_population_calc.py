
import numpy as np
from astropy import units
from astropy import constants

from prepare_data import *


def get_collision_rate_coefficient_excition(i,j, T, Eij, Yij, wi):
	'''
	collision rate per unit volume N_e*N_i*q_{ij} where N_i is number density of ion in level i, N_e is total number of free electron per unit volume
	q_{ij} is a measure of the frequency with which the free electron induce the transition i--j in palasam temperature T

	From "On the analysis of collision strengths and rate coefficients" Burgess and Tully 1992: eq(20)

	q_{ij} = 2*pi^{1/2}*a_0*hbar*m_e^{-1}*(I_inf/(kT))^{1/2}*exp(-E_{ij}/(kT))*Y_{ij}/w_i
	where Y_{ij} is thermally averaged collision strength, w_i is statistical weight of level i, E_{ij} is transition enerfy from level i to level j

	2*pi^{1/2}*a_0*hbar*m_e^{-1} = 2.1716*10**{-8} cm^3 s^{-1}

	INPUTS:
		i:
		j:
		T: temperature in Kelvin, number without unit
		Eij: transition energy in eV, number without unit
	'''

	Iinf = 33.50 #*units.eV   #??
	#Iinf = 13.6
	k = constants.k_B
	T = T*units.K	
	kT_eV = (k*T).to('eV').value

	#Eij = Eij*units.eV

	qij = 2.1716e-8*(Iinf/(kT_eV))**(1.0/2)*np.exp(-Eij/(kT_eV))*Yij/wi

	return qij

	
def get_collision_rate_coefficient_downward(T, Eij, qij, wi,wj):
	'''
	
	q_{ji} = (wi/wj)*exp(Eij/(kT))*q_{ij}
	'''

	k = constants.k_B
	
	T = T*units.K	
	#Eij = Eij*units.eV
	kT_eV = (k*T).to('eV').value	

	qji = (wi/wj)*np.exp(Eij/kT_eV)*qij

	#print qji
	return qji


def get_Avalue_array():

	A = np.zeros((15,15))
	for iindex in range(15):
		for jindex in range(15):
			#print iindex,jindex
			ilevel = iindex+1
			jlevel = jindex+1
			if ilevel < jlevel:
				Aji = get_transition_probability(jlevel, ilevel)
				A[jindex, iindex] = np.round(Aji,3)
			else:
				continue

	return A


def get_Eij_array():
	'''
	transition energy
	'''
	E = np.zeros((15,15))
	for iindex in range(15):
		for jindex in range(15):
			#print iindex,jindex
			ilevel = iindex+1
			jlevel = jindex+1
			if ilevel == jlevel:
				E[iindex,jindex]= 0
			elif ilevel < jlevel:
				Eij = get_transition_energy(ilevel, jlevel, eunit='eV')
				E[iindex, jindex] = np.round(Eij,3)
			else:
				continue

	return E
	

def get_Yij_array(T):
	'''
	thermally averaged collision strength
	'''
	Y = np.zeros((15,15))
	for iindex in range(15):
		for jindex in range(15):
			#print iindex,jindex
			ilevel = iindex+1
			jlevel = jindex+1
			if ilevel == jlevel:
				Y[iindex,jindex]= 0
			elif ilevel < jlevel:
				Yij = get_Tave_collision_strength(ilevel, jlevel, T)
				Y[iindex, jindex] = np.round(Yij,4)
			else:
				continue

	return Y

	

def get_collision_rate_coefficient_array(T, ):
	'''
	There shoule be two collision rate coefficient array, excitation (upward transition) and deexcitation (downward transition).

		
	'''

	Earray = get_Eij_array()
	Yarray = get_Yij_array(T)
	
	qarray_ex  = np.zeros((15, 15))
	qarray_dex = np.zeros((15, 15))	

	for iindex in range(15):
		for jindex in range(15):

			#print iindex,jindex
			ilevel = iindex+1
			jlevel = jindex+1

			wi = get_statistical_weight(ilevel)
			wj = get_statistical_weight(jlevel)	

			if ilevel < jlevel:
				qij_ex = get_collision_rate_coefficient_excition(ilevel,jlevel, T, Earray[iindex,jindex], Yarray[iindex, jindex], wi)
				qji_dex = get_collision_rate_coefficient_downward(T, Earray[iindex, jindex], qij_ex, wi,wj)

				qarray_ex[iindex, jindex] = qij_ex
				qarray_dex[jindex, iindex] = qji_dex
			else:
				continue	

	return qarray_ex, qarray_dex





if __name__ == "__main__":


	Te = 1e4	
	Ne = 1e7

	A = get_Avalue_array()
	AT = A.transpose()
	Asumrow = np.sum(A, axis=1)
	Asr_diag = np.diag(Asumrow)
	#print "A:"
	#print A
	#print "AT:"
	#print AT
	#print "Asumrow:"
	#print Asumrow
	#print "Asr_diag:"
	#print Asr_diag
	
	E = get_Eij_array()
	Y = get_Yij_array(Te)
	#print E
	#print Y

	qex, qdex = get_collision_rate_coefficient_array(Te)
	#print qex*Ne
	#print qdex*Ne

	qdexT = qdex.transpose() 
	qdexsumrow = np.sum(qdex, axis=1)
	qdexsr_diag = np.diag(qdexsumrow)

	#print qdexT
	#print qdexsumrow
	#print qdexsr_diag

	qexT = qex.transpose() 
	qexsumcol = np.sum(qex, axis=0)	
	qexsr_diag = np.diag(qexsumcol)

	print qex
	print qexsumcol
	print qexsr_diag

	uarray = np.ones((15,15))
	
	b = np.array([1]*15)

	ap = AT- Asr_diag + Ne*qdexT - Ne*qdexsr_diag - Ne*qex + Ne*qexsr_diag


	#print "decay:"
	#print AT-Asr_diag

	#print "de-excitation:"
	#print Ne*qdexT - Ne*qdexsr_diag


	#print "excitation:"
	#print -Ne*qex+Ne*qexsr_diag

	a = AT- Asr_diag + Ne*qdexT - Ne*qdexsr_diag - Ne*qex + Ne*qexsr_diag + uarray
	

	out = np.linalg.solve(a,b)
	print out

	sv = np.sum(out)
	print sv



