


def wavelength_convert_air_vacuum(lambda1, verbose=1):
	'''
	lambda_air = lambda_vac / n; where n is tmospheric refractivity

	n = 1 + 6.4328e-5 + (2.94981e6)/(1.46e10-sigma**2) + (2.5540e4)/(4.1e9-sigma**2); (old) where sigma is wave number

	n = 1 + 8.34213e-5 + (2.406030e6)/(1.30e10-sigma**2) + (1.5997e4)/(3.89e9-sigma**2)

	'''
	
	sigma = 1.0/lambda1*1e8

	#n = 1 + 6.4328e-5 + (2.94981e6)/(1.46e10-sigma**2) + (2.5540e4)/(4.1e9-sigma**2)
	n = 1 + 8.34213e-5 + (2.406030e6)/(1.30e10-sigma**2) + (1.5997e4)/(3.89e9-sigma**2)

	lambda2 = lambda1/n

	if verbose:
		print "The atmospheric refractivity at %s angstrom is %s"%(lambda1, n)

	return lambda2

if __name__ == "__main__":

	import sys
	lambda1 = float(sys.argv[1])

	outlambda = wavelength_convert_air_vacuum(lambda1)		
	print outlambda
	
