import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import warnings
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from tqdm import trange

	
eps = 1E-12

# Root mean squared error normalized by the mean of the observed data, y.
def CVRMSE(y, f):
	yavg = sum(y)/len(y)
	RMSE = np.sqrt( sum((y-f)**2) / len(y) )
	return ( abs(RMSE/yavg) )

# R squared for observed y, and predicted f.
def r2score(y, f):
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


def load_hashtable_to_dataframe(filename):
	''' For output of data from hashtables used in inversion_cascade_*.f90
	From a sparse data set of coordinates and other dictionary
	variables such as exits, probability, visits, etc. with metdata
	create a data frame and a dictionary of the metadata. '''
	metadata = {}
	array_keys = {'Burst', 'lmbda', 'k', 'n', 'c', 'beta', 'x_mean', 'r_mean'}
	data_lines = []
	header = None

	with open(filename, 'r') as f:
		for line in f:
			line = line.strip()
			if line.startswith('#'):
				# Parse metadata
				if ':' in line:
					key, value = line.lstrip('#').split(':', 1)
					key = key.strip()
					value = value.strip()
					parts = value.split()
					
					if key in array_keys:
						# Try float first, fallback to int
						try:
							metadata[key] = np.array([float(x) for x in parts])
						except ValueError:
							metadata[key] = np.array([int(x) for x in parts])
					else:
						# Try scalar float or int, fallback to string
						try:
							metadata[key] = float(value)
						except ValueError:
							metadata[key] = value
			elif line == '':
				continue
			else:
				if header is None:
					# First non-comment, non-blank line is the header
					header = line.split()
				else:
					data_lines.append(line.split())

	# Convert data to a DataFrame
	df = pd.DataFrame(data_lines, columns=header)
	df = df.apply(pd.to_numeric)

	return metadata, df


def compute_rates_df(this, time):
	'''Given exits probability and time, add rate columns to a df'''
	# Extract exit_fields for calling whenever we need them.
	exit_fields = [col for col in this.columns if col.startswith('exit')]

	for i, exit_col in enumerate(exit_fields, start=1):
		rate_col = f'rate{i}'
		this[rate_col] = this[exit_col] / (this['probability'] * time)
	
	return(this, exit_fields)


def forest_dependence_test(df, target_col, n_permutations=100, n_training=10, n_repeats=30, random_state=42):
	from sklearn.ensemble import RandomForestRegressor
	from sklearn.inspection import permutation_importance
	from sklearn.utils import shuffle

	""" Assess feature dependence using permutation importance and a null distribution.
	
	Returns: A DataFrame with importance, z-score, and empirical p-value per feature.
	"""
	rng = np.random.RandomState(random_state)
	X = df.drop(columns=target_col)
	y = df[target_col].values
	
	actual_importances = np.zeros((n_training, X.shape[1]))
	
	for i in range(n_training):
		# Train on real data
		model = RandomForestRegressor(n_estimators=200, min_samples_leaf=5, random_state=random_state+i)
		model.fit(X,y)
		result = permutation_importance(model, X, y, n_repeats=n_repeats, random_state=random_state)
		actual_importances[i,:] = result.importances_mean
		print(result.importances_mean)
		
	mean_importance = actual_importances.mean(axis=0)
	print('mean:', mean_importance)
	
	# Get null distribution
	null_importances = np.zeros((n_permutations, X.shape[1]))
	
	for i in range(n_permutations):
		y_perm = shuffle(y, random_state=rng)
		model_perm = RandomForestRegressor(n_estimators=200, min_samples_leaf=5, random_state=random_state+i)
		model_perm.fit(X,y_perm)
		perm_result = permutation_importance(model_perm, X, y_perm, n_repeats=n_repeats, random_state=random_state+i)
		null_importances[i,:] = perm_result.importances_mean
		print(perm_result.importances_mean)
	
	# Calculate z-scores and empirical p-values
	null_mean = null_importances.mean(axis=0)
	null_std = null_importances.std(axis=0)
	# z_scores = (actual_importance - null_mean) / (null_std + 1e-9)
	# p_values = np.mean(null_importances >= actual_importance[np.newaxis, :], axis=0)
	z_scores = (mean_importance - null_mean) / (null_std + 1e-9)
	p_values = np.mean(null_importances >= mean_importance[np.newaxis, :], axis=0)
	
	return pd.DataFrame({
		'feature': X.columns,
		# 'importance': actual_importance,
		'importance': mean_importance,
		'z_score': z_scores,
		'p_value': p_values
		}).sort_values(by='p_value')


  # Create a unique sentinel
_sentinel = object()
def permutation_importance_test(train_df, features, targets, n_permutations, test_df=_sentinel, test_size=0.3, seed=42):
	''' Asses feature dependence by comparing the performance of a
	random forest model trained on a fraction of our data to the
	remaining real data, and shuffled data. Idea is that if we can
	shuffle our data and the model still performs well, the feature we
	shuffled can not have been important to the target rate.
	
	INPUT our data as df, with features we want to test and target
	rates we are interested in.
	
	RETURN a dataframe of our target and features with their baseline
	error, mean_permutation error, std_permutation error, z_score, and
	p_score.'''
	
	rng = np.random.default_rng(seed)
	stats_dict = {}
	error_dict = {}
	
	for target in targets:
		X = train_df[features]
		y = train_df[target].to_numpy()
		# Make a test/training split if we don't provide training data
		if (test_df is _sentinel):
			X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)
		else:
			X_train = X
			y_train = y
			X_test = test_df[features]
			y_test = test_df[target].to_numpy()
		# Make the model
		model = RandomForestRegressor(n_estimators=200, random_state=seed)
		model.fit(X_train, y_train)
		# Get baseline prediction errors
		# baseline_preds = model.predict(X_test)
		# baseline_error = mean_squared_error(y_test, baseline_preds)
		baseline_error = model.score(X_test, y_test)
		
		for feature in features:
			perm_errors=[]			
			for _ in trange(n_permutations, desc=f"{target}: permuting {feature}", leave=True):
				# Copy and permute feature of interest
				X_perm = X_test.copy()
				X_perm[feature] = rng.permutation(X_perm[feature].values)
				# Test model against shuffled test data
				# preds = model.predict(X_perm)
				# err = mean_squared_error(y_test, preds)
				err = model.score(X_perm, y_test)
				perm_errors.append(err)
				
			# Calculate statistics
			perm_errors = np.array(perm_errors)
			mean_perm = np.mean(perm_errors)
			std_perm = np.std(perm_errors)
			z_score = (baseline_error - mean_perm) / std_perm if std_perm > 0 else np.nan
			p_value = np.mean(perm_errors >= baseline_error)
			
			# Keep statistics
			key = (target, feature)
			stats_dict[key] = {
				'baseline_error': baseline_error,
				'mean_perm_error': mean_perm,
				'std_perm_error': std_perm,
				'coeff_var': std_perm/mean_perm,
				'z_score': z_score,
				'p_value': p_value
			}
			
			# Keep errors in dictionary for histogram plotting.
			error_dict[(target, feature)] = perm_errors
	return (stats_dict, error_dict)


def plot_feature_permutation_histogram(ax, stats_dict, error_dict, target, feature, bins=30):
	"""
	Plot permutation error histogram for one target-feature on the given ax.
	
	Parameters:
	- ax: matplotlib Axes object to plot on
	- stats_df: DataFrame with baseline errors (indexed by target)
	- error_dict: dict with keys (target, feature), values = permutation error arrays
	- target: string, target variable name
	- feature: string, feature name
	- bins: int or array, histogram bins
	"""
	errors = error_dict.get((target, feature), None)
	ax.hist(errors, bins=bins, density=True)

	# Plot baseline error from stats_df
	baseline = stats_dict[target, feature]['baseline_error']
	ax.axvline(baseline, color='k', linestyle='--', linewidth=2, label='Baseline')

	ax.set_title(f"Permutation Errors for {feature}")
	# ax.set_xlabel("Error")
	# ax.set_ylabel("Density")

