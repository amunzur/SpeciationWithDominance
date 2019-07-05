#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_file(n, N, alpha, u, sigma, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f' %(n, N, alpha, u, sigma)
	outfile_A = open("%s/PlotBurn_%s.csv" %(data_dir, sim_id), "w")
	return outfile_A

def write_data_to_output(fileHandles, data):
	"""
	This function writes to the corresponding output file.
	"""
	writer = csv.writer(fileHandles)
	writer.writerow(data)

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""
	fileHandles.close()

def save_arrays(n, N, alpha, u, sigma, rep, data_dir, mut, pop):
	"""
	Save numpy arrays of mutations and their frequencies.
	"""
	sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f_rep%d' %(n, N, alpha, u, sigma, rep)
	
	filename = "%s/Muts_%s.npy" %(data_dir, sim_id)
	np.save(filename, mut[1:]) #save mutations
	
	filename = "%s/Freqs_%s.npy" %(data_dir, sim_id)
	np.save(filename, np.sum(pop[:,1:], axis=0)/len(pop)) #save frequencies

def fitness(phenos, theta, sigma):
	"""
	This function determines relative fitness
	"""
	dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
	w = np.exp(-0.5 * sigma * dist**2) #fitness
	return w

def recomb(surv):
	"""
	This function creates offspring through pairing of parents (haploid) and recombination (i.e, meiosis)
	"""
	pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	off = np.append(off_1, off_2, axis=0) #each product of meiosis
	return off

def mutate(off, u, mean_mut, cov_mut, mut):
	"""
	This function creates mutations and updates population
	"""
	rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
	nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
	whomuts = np.where(rand3 < u) #indices of mutants
	newmuts = np.random.multivariate_normal(mean_mut, cov_mut, nmuts)
	# newmuts = np.random.normal(0, alpha, (nmuts, m)) #phenotypic effect of new mutations
	pop = np.append(off, np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
	mut = np.append(mut, newmuts, axis=0) #append effect of new mutations to mutation list
	return [pop, mut]

def remove_lost_muts(remove_lost, pop, mut):
	"""
	This function removes lost mutations
	"""
	if remove_lost:
		keep = pop.any(axis=0)
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

def remove_fixed_muts(remove_fixed, pop, mut):
	"""
	This function removes fixed mutations
	"""
	if remove_fixed:
		keep = np.concatenate((np.array([True]), np.any(pop[:,1:]-1, axis=0))) #note this is a little more complicated because we want to keep the first, base, mutation despite it being fixed
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

def histogram_files(pop, mut, theta, mean_mut, cov_mut, n, N, alpha, u, sigma, rep, data_dir):
	"""
	Save csv of mutation sizes (in SGV and de novo) for plotting histograms
	"""
	if make_histogram_files:

		#segregating mutations
		keep = np.any(pop-1, axis=0)

		sgv_dist = np.linalg.norm(mut[keep] - theta, axis=1) #phenotypic distance from optimum for each individual mutation in sgv
		sgv_freq = np.sum(pop[:, keep], axis=0) #number of copies of each mutation
		newmuts = np.random.multivariate_normal(mean_mut, cov_mut, len(sgv_dist*10)) #phenotypic effect of new mutations (make 10 times number as in sgv to get clean distribution)
		dist_denovo = np.linalg.norm(newmuts - theta, axis=1) #phenotypic distance from optimum for each individual de novo mutation

		sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f_rep%d' %(n, N, alpha, u, sigma, rep) #sim info
		filename_sgv = "%s/sgv_muts_%s.csv" %(data_dir, sim_id) #filename for sgv mutations
		filename_denovo = "%s/denovo_muts_%s.csv" %(data_dir, sim_id) #filename for de novo mutations

		#write sgv csv
		with open(filename_sgv, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerow(sgv_dist)
			writer.writerow(sgv_freq)

	    #write de novo csv
		with open(filename_denovo, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerow(dist_denovo)

		#distribution of trait values among individuals
		# z = np.dot(pop, mut) #phenotypes		
		# filename_z = "%s/z_%s.csv" %(data_dir, sim_id) #filename
		# # fileHandle_B(filename_z, z)
		# with open(filename_z, 'w') as csvfile:
		# 	writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		# 	writer.writerow(z)

######################################################################
##PARAMETERS##
######################################################################

ns = [2] #phenotypic dimensions (positive integer >=1)
Ns = [10000] #number of haploid individuals (positive integer >=1)
u = 0.001 #mutation probability per generation per genome (0<=u<=1)
sigmas = [0.01] #strength of selection (positive real number)
alpha = 0.1 #mutational sd in each trait dimension (positive real number)

maxgen = 25000 #total number of generations population adapts for (positive integer)
gen_rec = 500 #print every this many generations (positve integer <=maxgen)

remove_lost = True #If true, remove mutations that are lost
remove_fixed = False #If true, remove mutations that are fixed

make_histogram_files = True #if true ouput mutation sizes for plotting 

reps = 5 #number of replicates (positive integer)

data_dir = 'data/burnin' #where to save data

######################################################################
##MAIN SIMULATION CODE##
######################################################################

def main():

	#loop over population size
	i_N = 0
	while i_N < len(Ns):
		N = Ns[i_N]

		#loop over selection strength
		i_sigma = 0
		while i_sigma < len(sigmas):
			sigma = sigmas[i_sigma]

			#loop over dimensions
			i = 0
			while i < len(ns):
				n = ns[i]
				
				#set the n-dependent parameters
				mean_mut = [0] * n #mean mutational effect in each dimension (real^n)
				cov_mut = np.identity(n) * alpha**2 #mutational covariance matrix (real^nxn)
				theta = np.array([0] * n) # optimum phenotype (n real numbers)

				#open files for saving data
				fileHandles = open_output_file(n, N, alpha, u, sigma, data_dir)

				#loop over replicates
				rep = 1
				while rep < reps + 1:

					#initialize the population
					pop = np.array([[1]] * N) #start all individuals with one "mutation"
					mut = np.array([[0] * n]) #but let the "mutation" do nothing (ie put all phenotypes at origin)

					#intitialize generation counter
					gen = 0

					#run until maxgen
					while gen < maxgen + 1:

						# genotype to phenotype
						z = np.dot(pop, mut) #sum mutations held by each individual (ie additivity in phenotype)

						# phenotype to fitness
						w = fitness(z, theta, sigma)

						# wright-fisher (multinomial) sampling
						parents = np.random.multinomial(N, w/sum(w)) #number of times each parent chosen
						parent_phenos = np.dot(pop[parents], mut)	
						psegvar = np.mean(np.var(parent_phenos, axis = 0)) # segregation variance (mean of individual trait variances)
						off = np.repeat(pop, parents, axis=0) #offspring genotypes

						# mating and recombination
						off = recomb(off)

						# mutation and population update
						[pop, mut] = mutate(off, u, mean_mut, cov_mut, mut)

						# remove lost mutations if doing so
						[pop, mut] = remove_lost_muts(remove_lost, pop, mut)

						#print update every gen_rec generations	
						if gen % gen_rec == 0:
							
							# print(np.mean(pop[:,1:],axis=0))
							print('N = %d, sigma = %.2f, n = %d   rep=%d   gen = %d   num_seg=%d   mean_freq=%.3f   mean_dist=%.3f   psegvar=%.3f' %(N, sigma, n, rep, gen, len(mut), np.mean(np.mean(pop[:,1:],axis=0)), np.mean(np.linalg.norm(mut - theta, axis=1)), psegvar))

							#save for plotting approach to MS balance (number of segregating mutations, avg frequency)
							seg = np.any(pop-1, axis=0) #segregating mutations only
							mut_seg = mut[seg] #which mutations segregating
							pop_seg = pop[:, seg] #only look at those sites which are segregating
							write_data_to_output(fileHandles, [rep, gen, len(mut_seg), np.mean(np.mean(pop_seg,axis=0)), np.mean(np.linalg.norm(mut_seg - theta, axis=1)), psegvar]) #save replicate, generation, number of segregating mutations, avg frequency, and average departure from optimum

						# go to next generation
						gen += 1

					# remove fixed mutations if doing so
					[pop, mut] = remove_fixed_muts(remove_fixed, pop, mut)

					#save mutation and frequency data
					save_arrays(n, N, alpha, u, sigma, rep, data_dir, mut, pop)

					#save mutation sizes for plotting histogram
					histogram_files(pop, mut, theta, mean_mut, cov_mut, n, N, alpha, u, sigma, rep, data_dir)

					# go to next rep
					rep += 1

				#clean up
				close_output_files(fileHandles)

				#go to next n value
				i += 1

			#go to next sigma value
			i_sigma += 1

		#go to next N value
		i_N += 1

######################################################################
##RUN SIMULATIONS##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
