#! /anaconda/bin/python

def get_superscript_number(n):
	'''
	see https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts

	INPUTS:
		n: 
	'''
	
	supdict = {'0':u'\u2070', '1':u'\u00B9', '2':u'\u00B2', '3': u'\u00B3', '4':u'\u2074', '5':u'\u2075', '6':u'\u2076', '7':u'\u2077', '8':u'\u2078', '9':u'\u2079',}
	out = supdict[str(n)]

	return out


def get_subscript_number(n):
	'''
	INPUTS:
		n: 
	'''
	subdict = {'0':u'\u2080', '1':u'\u2081', '2':u'\u2082', '3': u'\u2083', '4':u'\u2084', '5':u'\u2085', '6':u'\u2086', '7':u'\u2087', '8':u'\u2088', '9':u'\u2089',}
	out = subdict[str(n)]
	return out

def electron_configuration(n,l,ne):
	'''
	INPUTS:
		n: principal quant number
		l: orbital quant number
		ne: electron number
	'''
	
	supne = get_superscript_number(ne)
	out = u'%s%s%s'%(n,l,supne)
	
	return out

def spectrum_term(multiplicity, L, J_numerator, J_dominator):
	'''
	INPUTS:
		multiplicity: 2S+1,  
		L: total orbital angular momentum
		J_numerator: 
		J_dominator:
	'''

	supm = get_superscript_number(multiplicity)	
	subJn = get_subscript_number(J_numerator)
	subJd = get_subscript_number(J_dominator)
	slash = u'\u2044'

	out = u'%s%s%s%s%s'%(supm, L, subJn, slash, subJd)
	#out = supm + u'%s'%L + subJn + slash + subJd

	return out

if __name__ == "__main__":
	
	example1 = spectrum_term(4,'F',9,2)
	print example1

	example2 = electron_configuration(3,'d',7)
	print example2
