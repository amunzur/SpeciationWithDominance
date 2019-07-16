#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
# import matplotlib.pyplot as plt
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, N, alpha, u, sigma, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_N%d_alpha%.4f_u%.4f_sigma%.4f' %(n, N, alpha, u, sigma)
	outfile_A = open("%s/Fig3_4_%s.csv" %(data_dir, sim_id), "w") #summary data
	outfile_B = open("%s/Fig3_4_phenotypes_%s.csv" %(data_dir, sim_id),"wb") #hybrid phenotypes
	return [outfile_A, outfile_B]

def write_data_to_output(fileHandles, data):
	"""
	This function writes a (time, data) pair to the
	corresponding output file. We write densities
	not abundances.
	"""
	writer = csv.writer(fileHandles)
	writer.writerow(data)

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""
	fileHandles.close()

def fitness(phenos, theta, sigma):
	"""
	This function determines relative fitness
	"""
	dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
	w = np.exp(-0.5 * sigma * dist**2) #fitness
	return w

def shuffle_along_axis(a, axis):
	idx = np.random.rand(*a.shape).argsort(axis = axis)
	return np.take_along_axis(a,idx,axis = axis) 


def recomb(surv):
	"""
	This function creates offspring through pairing of parents (diploid) and recombination (i.e, meiosis)
	"""

	# crossing over within diploid parents to make gametes, input is 1000 rows, output also 1000 rows; n_loci columns
	surv_stacked = np.stack((off1, off2), axis = 1).reshape(N_adapts * 2, 4) # this places the related rows together, 1st row of each together, then 2nd, then 3rd... - 4 loci 
	
	#recombination 
	c = np.array([0,0,0,0,]).reshape([1, 4]) # create a random array to save the results of each loop

	x = 0
	while x < 21:
		b = shuffle_along_axis(surv_stacked[x:(x+2)], axis = 0) #shuffle along columns in groups of 2. each group of 2 represents chrom1 an chrom2 of each individual 
		surv_stacked = np.concatenate((c, b), axis = 0) #update what surv_stacked is after each loop of recombination 
		x+=2

	surv_stacked = surv_stacked[1:(N_adapts + 1)] #remove the c from the top of the array, update surv_stacked accordingly 

	surv_chrom1 = surv_stacked[::2] #this selects every odd row - chrom1 of N_adapts number of individuals 
	surv_chrom2 = surv_stacked[1::2] #this selects every even row - chrom 2 
	
	# pick random pairs of parents 
	pairs = np.resize(np.random.choice(len(surv_stacked), size=len(surv_stacked), replace=False), (int(len(surv_stacked)/2), 2))

	#determine the gamates of parents 
	parent1_meiosis1 = surv_chrom1[pairs[:, 0]] #pick the chrom1 of the parents (pick the related rows, 0th element of pairs )
	parent1_meiosis2 = surv_chrom2[pairs[:, 0]]

	parent2_meiosis1 = surv_chrom1[pairs[:, 1]]
	parent2_meiosis2 = surv_chrom2[pairs[:, 1]]

	# shuffle the gametes into new diploid individuals


	# output off array (matrix)


	# RANDOMLY PAIR THE PARENTAL CHROMOSOMES TO MAKE THE OFFSRPING 

	l = int(np.random.randint(low = 0, high = 4)) # randomly picks 0, 1, 2, 3 to randomly generate offsptring out of 4 possibilities 

	if l == 0:
		off = np.stack((parent1_meiosis1, parent2_meiosis1), axis = 1)
	
	elif l == 1:
		off = np.stack((parent1_meiosis1, parent2_meiosis2), axis = 1)
	
	elif l == 2:
		off = np.stack((parent1_meiosis2, parent2_meiosis1), axis = 1)
	
	else:	
		off = np.stack((parent1_meiosis2, parent2_meiosis2), axis = 1)

	return off # we made the offspring, now they will act as parents and we will create mutations in them. 


	# ORIGINAL CODE: 
	# pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	# rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	# rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	# off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	# off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	# off = np.append(off_1, off_2, axis=0) #each product of meiosis, diploid offspring
	# return off


def mutate(off, u, alpha, n, mut):
	"""
	This function creates mutations and updates population
	"""
	


	

	# h == 0 #looping over all the loci in off 
	# while h < len(off): 
	# 	h +=1
	# 	rand3 = np.random.uniform(size = len(off))
	# 	off[off > rand3] = 1
		


	# 	if rand < np.random.uniform(size = 1) is True: #if the loci is likely to mutate







	#ORIGINAL CODE
	rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring [or loci??]. (creates number of offxrandom numbers as between 0 and 1)
	nmuts = sum(rand3 < u_adapt) # mutate if random number is below mutation rate; returns number of new mutations - u is the mutation probability per generation per genome - ancestor
	whomuts = np.where(rand3 < u_adapt) #indices of mutants, meaning where they are 
	newmuts = np.random.normal(0, alpha, size = (nmuts, n)) #phenotypic effect of new mutations. 0=mean, alpha=sd rows:nmuts coloumns=n 
	pop = np.append(off, np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
	mut = np.append(mut, newmuts, axis=0) #append effect of new mutations to mutation list
	return [pop, mut]

def remove_muts(remove, remove_lost, pop, mut, mutfound):
	"""
	This function creates mutations and updates population
	"""
	if remove_lost:
		if remove == 'any':
			keep = pop.any(axis=0)
			mut = mut[keep]
			pop = pop[:, keep]
		elif remove == 'derived':
			segregating = pop.any(axis=0)
			ancestral = np.array(range(len(mut))) < len(mutfound)
			keep = np.add(segregating, ancestral)
			mut = mut[keep]
			pop = pop[:, keep]
	return [pop, mut]

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

nreps = 10 #number of replicates for each set of parameters
ns = [2] #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

N_adapts = [1000] #number of diploid individuals (positive integer)
alpha_adapt = 0.1 #mutational sd (positive real number)
u_adapt = 0.001 #mutation probability per generation per genome (0<u<1)
sigma_adapts = [10] #selection strengths

opt_dist = 1 #distance to optima

n_angles = 3 #number of angles between optima to simulate (including 0 and 180) (>=2)

maxgen = 50 #total number of generations populations adapt for

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##PARAMETERS FOR HYBRIDS##
######################################################################

nHybrids = 100 #number of hybrids to make at end of each replicate

######################################################################
##FUNCTION FOR POPULATIONS TO ADAPT##
######################################################################

def main():

	#loop over population size
	i_N = 0
	while i_N < len(N_adapts):
		N_adapt = N_adapts[i_N]

		#loop over selection strength
		i_sigma = 0
		while i_sigma < len(sigma_adapts):
			sigma_adapt = sigma_adapts[i_sigma]

			#loop over dimensions
			l = 0
			while l < len(ns):
				n = ns[l]

				#set n-dependent parameters
				theta1 = np.append(opt_dist,[0]*(n-1)) #set one optima to be fixed
				angles = [math.pi*x/(n_angles-1) for x in range(n_angles)] #angles to use (in radians)
				if n == 2:
					theta2_list = np.array([[opt_dist*math.cos(x), opt_dist*math.sin(x)] for x in angles]) #optima to use
				elif n > 2:
					theta2_list = np.array([np.append([opt_dist*math.cos(x), opt_dist*math.sin(x)], [0]*(n-2)) for x in angles]) #optima to use

				# open output files
				[fileHandle_A, fileHandle_B] = open_output_files(n, N_adapt, alpha_adapt, u_adapt, sigma_adapt, data_dir)

				#loop over optima
				j = 0
				while j < len(theta2_list):
					
					#set optima
					theta2 = theta2_list[j]

					# #set up plot of hybrid load versus number of ancestral mutations (n_muts)
					# plt.axis([0, max(n_mut_list)+1, 0, 0.1])
					# plt.ylabel('hybrid load at generation %d (mean $\pm$ SD of %d replicates)' %(maxgen,nreps))
					# plt.xlabel('number of ancestral mutations')
					# plt.ion()

					#loop over all replicates
					rep = 0
					while rep < nreps:

						#found identical populations
						popfound1 = np.array([[1]] * N_adapt) #this creates an array of 1 coloumn of ones and N_adapt times rows. 
						popfound2 = np.array([[1]] * N_adapt) 
						popfound = np.column_stack((popfound1, popfound2)) # create a diploid genotype. #of columns= #of chromosomes 

						mutfound1 = np.array([[0]] * n) #similar to above. filled with zeroes. number of coloumns: n. rows: 1 
						mutfound2 = np.array([[0]] * n)
						mutfound = np.column_stack((mutfound1, mutfound2))

						[pop1, mut1] = [popfound, mutfound] # this creates pop and mut arrays for both parents. they are the same because we start from the same point. 
						[pop2, mut2] = [popfound, mutfound] # mut1 = how farther you go from the origin due to mutations in pop1. same for mut2

						pop1_chrom1 = pop1 # genotype of 1st chromosome of pop1
						pop1_chrom2 = pop1 # genotype of 2nd chromosome of pop1

						pop2_chrom1 = pop2 # 1st chromosome of pop2
						pop2_chrom2 = pop2 # 2nd chromosome of pop2

						# intitialize generation counter
						gen = 0

						# run until maxgen
						while gen < maxgen + 1:

							# genotype to phenotype (haploid):
							# phenos1 = np.dot(pop1, mut1) #sum mutations held by each individual (phenotype of the mutations)
							# phenos2 = np.dot(pop2, mut2) #sum mutations held by each individual 

							# genotype to phenotype (diploid):
							pop1_overall = ((pop1_chrom1 + pop1_chrom2) / 2) # two chromosomes of pop1 averaged
							pop2_overall = ((pop2_chrom1 + pop1_chrom2) / 2) # two chromosomes of pop2 averaged

							phenos1 = np.dot(pop1_overall, mut1)
							phenos2 = np.dot(pop2_overall, mut2)

							# phenotype to fitness
							w1 = fitness(phenos1, theta1, sigma_adapt)
							w2 = fitness(phenos2, theta2, sigma_adapt)

							# wright-fisher (multinomial) sampling
							parents1 = np.random.multinomial(N_adapt, w1/sum(w1)) # number of times each parent chosen, drawing samples from a multinomial ditribution
							# N_adapt = number of experiments, w1/sum(w1 = probability of parent1 being chosen. if you are more fit, you are chosen more often. 
							
							off1_chrom1 = np.repeat(pop1_overall, parents1, axis=0) 
							off1_chrom2 = np.repeat(pop1_overall, parents1, axis=0)
							off1 = [off1_chrom1, off1_chrom2] 
							
							parents2 = np.random.multinomial(N_adapt, w2/sum(w2)) # number of times each parent chosen
							off2_chrom1 = np.repeat(pop2_overall, parents2, axis=0) # offspring genotypes of pop2
							off2_chrom2 = np.repeat(pop2_overall, parents2, axis=0)
							off2 = [off2_chrom1, off2_chrom2]

							# mating and recombination
							off1 = recomb(off1)
							off2 = recomb(off2)


							# mutation and population update
							[pop1_overall, mut1] = mutate(off1, u_adapt, alpha_adapt, n, mut1)
							[pop2_overall, mut2] = mutate(off2, u_adapt, alpha_adapt, n, mut2)

							# remove lost mutations (all zero columns in pop)
							[pop1_overall, mut1] = remove_muts(remove, remove_lost, pop1_overall, mut1, mutfound)
							[pop2_overall, mut2] = remove_muts(remove, remove_lost, pop2_overall, mut2, mutfound)

							# go to next generation
							gen += 1

						#parent fitness and load (use parent 1, but could be either)	
						parents = np.random.randint(len(pop1), size = nHybrids)
						parent_phenos = np.dot(pop1_overall[parents],mut1)	
						mean_parent_pheno = np.mean(parent_phenos, axis=0)
						parent_fitnesses = fitness(parent_phenos, mean_parent_pheno, sigma_adapt) #parent fitnesses
						pfit = np.mean(parent_fitnesses) #mean parent fitness

						#make variables to hold offspring phenotypes
						offphenos = dict()
						offpheno = []

						#make each of nHybrids hybrids
						for k in range(nHybrids):
							# choose random parents
							randpar1 = pop1_overall[np.random.choice(len(pop1_overall))] 
							randpar2 = pop2_overall[np.random.choice(len(pop2_overall))]
							# get random parent phenotypes
							phenpar1 = np.dot(randpar1, mut1)
							phenpar2 = np.dot(randpar2, mut2)
							# get mutations held by random parents
							mutpar1 = mut1 * randpar1[:, None]
							mutpar2 = mut2 * randpar2[:, None]
							setA = set(tuple(x) for x in mutpar1)
							setB = set(tuple(x) for x in mutpar2)
							# find mutations shared by two parents (all in offspring)
							sharedmuts = np.array([x for x in setA & setB])
							if len(sharedmuts) < 1:
								sharedmuts = np.array([[0] * n]) #give something in case empty
							# find mutations not shared by two parents
							unsharedmuts = np.array([x for x in setA ^ setB])
							# which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
							randmuts = np.random.randint(2, size = (len(unsharedmuts)))	
							unsharedoffmuts = unsharedmuts * randmuts[:, None]
							if len(unsharedoffmuts) < 1:
							    unsharedoffmuts = np.array([[0] * n]) #give something in case empty
							# offspring phenotype is collection of shared and random unshared mutations
							offpheno.append(sum(np.append(sharedmuts, unsharedoffmuts, axis = 0)))

						#mean hybrid fitness
						hybrid_fitnesses = np.maximum(fitness(offpheno, theta1, sigma_adapt), fitness(offpheno, theta2, sigma_adapt)) #max fitness of hybrid, ie. gets fitness in parent enviro it is closest to
						hyfit = np.mean(hybrid_fitnesses) #mean fitness
						rhyfit = hyfit/pfit #relative fitness
						max_hyfit = np.percentile(hybrid_fitnesses, 95) # max fitness over all hybrids (90th percentile ie top 10 per cent)
						rel_max_hyfit = max_hyfit/pfit #max hybrid fitness relative to parental mean

						#lag load
						meanpheno = np.mean(offpheno, axis=0) #mean hybrid phenotype
						fitmeanpheno = np.mean(np.maximum(fitness(np.array([meanpheno]), theta1, sigma_adapt), fitness(np.array([meanpheno]), theta2, sigma_adapt)))/pfit #fitness of mean hybrid phenotype
						# lagload = -np.log(fitmeanpheno) #lag load
						Emeanpheno = np.mean(np.array([theta1, theta2]), axis=0)
						Efitmeanpheno = np.mean(fitness(np.array([Emeanpheno]), theta1, sigma_adapt))/pfit

						#save hybrid phenotype data (to make hybrid clouds)
						pheno_data = np.column_stack( [np.array([i+1 for i in range(len(offpheno))]), np.array([rep+1]*len(offpheno)), np.array([round(angles[j]*180/math.pi,2)]*len(offpheno)), np.array([opt_dist * (2*(1-math.cos(angles[j])))**(0.5)]*len(offpheno)), offpheno]) #make hybrid phenotype data (with a bunch of identifiers in front of each phenotype)
						np.savetxt(fileHandle_B, pheno_data, fmt='%.3f', delimiter=',') #save

						#save summary data
						write_data_to_output(fileHandle_A, [round(angles[j]*180/math.pi,2), rep+1, opt_dist * (2*(1-math.cos(angles[j])))**(0.5), rhyfit])

						#print an update
						print('N = %d, sigma = %.2f, n=%d, angle=%r, rep=%d' %(N_adapt, sigma_adapt, n, round(angles[j]*180/math.pi,2), rep+1)) 

						# go to next rep
						rep += 1

						# #plot mean and SD hybrid load over all replicates for this n_muts value
						# plt.errorbar(n_mut_list[i], np.mean(hyloads), yerr=np.var(hyloads)**0.5, fmt='o', color='k')
						# plt.pause(0.01)

					#go to next optimum
					j += 1

				# cleanup
				close_output_files(fileHandle_A)
				close_output_files(fileHandle_B)

				#next dimension
				l += 1

			#go to next sigma value
			i_sigma += 1

		#go to next N value
		i_N += 1

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    

start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))