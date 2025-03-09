import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import warnings

eps = 1E-12

def CVRMSE(y, f):
# Root mean squared error normalized by the mean of the observed data, y.
	yavg = sum(y)/len(y)
	RMSE = np.sqrt(sum((y-f)**2)/len(y))
	return ( abs(RMSE/yavg) )

def r2score(y, f):
# R squared for observed y, and predicted f.	
	yavg = sum(y)/len(y)
	ssres = sum((y-f)**2)
	sstot = sum((y-yavg)**2)
	return(1 - ssres/sstot)
	
def relative_SSE(y, f):
	mean = sum(y*f)/len(y)
	return( SSE(y, f) / mean )

# Sum of Squared Errors
def SSE(y, f):
	return (sum((y-f)**2))
	
# Sum of Errors
def SE(y,f):
	return (sum(y-f))

def sigmoid(x):
	return(1/(1 + np.exp(-x)))
	
def hill_pos(x,k,n):
	return(x**n / (x**n + k**n))


def reorder_vectors(vector1, vector2):
    # Zip two vectors together, sort based on vector1, and unzip the result
    sorted_pairs = sorted(zip(vector1, vector2))  # Sort based on the first vector

    # Unzip into two separate sorted lists
    sorted_vector1, sorted_vector2 = zip(*sorted_pairs)

    # Convert the tuples back to lists if needed
    return list(sorted_vector1), list(sorted_vector2)

# Read metadata for older files
def read_metadata(fname):
	metadata = {}
	# Open the file and read line by line
	with open(fname, "r") as file:
		for line in file:
			# Check if the line contains parameter metadata
			if line.startswith(" # Parameter Metadata"):
				# Start reading parameter metadata
				for line in file:
					if line.startswith(" # Data"):
						# Stop reading parameter metadata if "Data" section is encountered
						break
					# Extract parameter name and value
					if ":" in line:
						key, value = line.strip().split(":", 1)
						key = key.lstrip("#").strip()
						value = float(value)
						metadata[key.strip()] = value
	return(metadata)


def save_custom_pdf_metadata(fig_obj, target_path, metadata):
	## Enter the context manager for the pdf object
	with PdfPages(target_path) as pdf_obj:
		## Get metadata InfoDict object
		infodict = pdf_obj.infodict()
		## Suppress warnings, but ONLY for the metadata addition step
		with warnings.catch_warnings():
			## Iterate over metadata elements and add them to infodict
			for key, value in metadata.items():
				infodict[key] = value
		## Save the current figure
		pdf_obj.savefig(fig_obj, bbox_inches='tight')
		## Watch out for:
		## - more possible kwargs in the call to savefig
		## and/or
		## - needs to be filtered by format at a higher level, as this only works for pdf!
	## 
	return fig_obj


# Parse newer files. The files now include the correlation data, several
# means and etas, and the probability matrix.
def parse_file(filename):
	metadata = {}
	data = []
	probability_matrix = []
	headers = []
	
	with open(filename, 'r') as file:
		section = None
		for line in file:
			line = line.strip()
			
			# Skip vertical whitespace
			if not line:
				continue
			
			# Check which section we are in
			if line.startswith('#'):
				if line.startswith('# Parameter Metadata'):
					section = 'metadata'
					continue
				elif line.startswith('# Data'):
					section = 'data'
					continue
				elif line.startswith('# Probability Matrix'):
					section = 'probability_matrix'
					continue
				
			# Process lines based on the current section	
			if section =='metadata' and ':' in line:
				key, value = line.split(':', 1)
				key = key.lstrip('#').strip()
				value = float(value)
				metadata[key] = value
			elif section == 'data':
				if line.startswith('#') and not headers:
					headers = line.lstrip('#').strip().split()
				else:
					data.append([float(x) for x in line.split()])
			elif section == 'probability_matrix':
				probability_matrix.append([float(x) for x in line.split()])
	
	# Convert data and probability_matrix to numpy arrays
	data = np.array(data)
	probability_matrix = np.array(probability_matrix)
	return(metadata, data, probability_matrix)


def read_datafile(fname):
	metadata, data, probability = parse_file(fname)
		
	tsim = data[:,0]
	Amm = data[:,1]
	Apm = data[:,2]
	Amp = data[:,3]
	App = data[:,4]
	Arpbp = data[:,5]
	Arpdp = data[:,6]
	dt = tsim[1]-tsim[0]
	
	# Get mean to whatever power.
	abund = np.arange(0, np.max(probability.shape), step=1)
	mprob = np.sum(probability, axis=0)
	pprob = np.sum(probability, axis=1)
	mlavg = np.dot(mprob, np.power(abund[:len(mprob)], metadata['l1']))
	plavg = np.dot(pprob, np.power(abund[:len(pprob)], metadata['l2']))
	
	# Protein lifetime
	tau_p = metadata['pavg'] / metadata['Rpdavg']
		
	return (metadata, data, probability,
			tsim, Amm, Apm, Amp, App, Arpbp, Arpdp, dt,
			plavg, tau_p
			)


# Theoretical etas and means for the simple mrna-protein system
def linear_thry(t, params):
	beta = [params['beta_m'], params['beta_p']]
	lmbda = params['lmbda']
	alpha = params['alpha']

	mean = np.zeros(2)
	mean[0] = lmbda/beta[0]
	mean[1] = mean[0]*alpha/beta[1]

	eta = np.zeros((2,2))
	eta[0,0] = 1./mean[0]
	eta[0,1] = eta[0,0] * beta[1]/sum(beta)
	eta[1,0] = eta[0,1]
	eta[1,1] = eta[0,1] + 1./mean[1]
	
	A = np.zeros((4,len(t)))
	
	# Amm
	A[0,:] = eta[0,0]*mean[0]*mean[0] * np.exp(-t*beta[0])
	
	# Apm
	if (beta[0] != beta[1]):
		A[1,:] = (-np.exp(-t*beta[1])*eta[0,1]*mean[0]*mean[1] + 
			(-np.exp(-t*beta[1]) + np.exp(-t*beta[0]))*lmbda/beta[0]/beta[1] * eta[0,0]*mean[0]*mean[0] / 
			(1/beta[0] - 1/beta[1]))
	else:
		A[1,:] = (np.exp(-t*beta[1])*eta[0,1]*mean[0]*mean[1] + 
			t*np.exp(-t*beta[1])*lmbda/beta[0]*beta[1]*eta[0,0]*mean[0]*mean[0])
	
	# Amp
	A[2,:] = eta[0,1]*mean[0]*mean[1] * np.exp(-t*beta[0])
		
	# App
	if (beta[0] != beta[1]):
		A[3,:] = (-np.exp(-t*beta[1])*eta[1,1]*mean[1]*mean[1] + 
			(-np.exp(-t*beta[0]) + np.exp(-t*beta[1]))*lmbda/beta[0]/beta[1] * eta[0,1]*mean[0]*mean[1] / 
			(1/beta[0] - 1/beta[1]))
	else:
		A[3,:] = (np.exp(-t*beta[1])*eta[1,1]*mean[1]*mean[1] + 
			t*np.exp(-t*beta[1])*lmbda/beta[0]*beta[1]*eta[0,1]*mean[0]*mean[1])
			
	return(mean, eta, A)
	
	
def funcL1(l2, etapp, mavg, pavg, tp, tm):
	f = np.sqrt((etapp-1/(l2*pavg)) * l2*mavg * ((tp+l2*tm)/tm))
	return (f)

