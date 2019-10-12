#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
# import matplotlib.pyplot as plt
import csv
import math
import pandas as pd
import os
from tqdm import tqdm
import progressbar as pb


######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, N, alpha, u, sigma, opt_dist, dom):
	"""
	This function opens the output files and returns file
	handles to each.
	"""

	data_dir = '/Users/student/Desktop/asli/summer-2019/data_done'

	sim_id = '_n%s_N%s_alpha%s_u%s_sigma%s_opt_dist%s_dom%s' %(n, N, alpha, u, sigma, opt_dist, dom)
	outfile_A = open(os.path.join(data_dir, "mutsummary%s.csv" %(sim_id)), "w") #summary data
	outfile_B = open(os.path.join(data_dir, "phenotypes%s.csv" %(sim_id)), "w") #hybrid phenotypes at the end 
	outfile_C = open(os.path.join(data_dir, "adaptation%s.csv" %(sim_id)), "w")
	outfile_D = open(os.path.join(data_dir, "fitness%s.csv" %(sim_id)), "w")#save the fitness of the parents and the hybrids in the same file 
	outfile_E = open(os.path.join(data_dir, "fitMean%s.csv" %(sim_id)), "w")

	return [outfile_A, outfile_B, outfile_C, outfile_D, outfile_E]

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

def crossover(chromosome):
	
	c = np.random.randint(1, size = np.size(chromosome, 1)).reshape(1, np.size(chromosome, 1)) # create a random array to save the results of each loop
	
	x = 0
	while x < (N_adapts[0] * 2 + 1): 
		
		b = shuffle_along_axis(chromosome[x:(x+2)], axis = 0) #shuffle along columns in groups of 2. each group of 2 represents chrom1 an chrom2 of each individual 
		c = np.concatenate((c, b), axis = 0) #update what F1_after_recomb2 is after each loop of recombination 
		x+=2

	crossedover = c[1:(N_adapts[0] * 2 + 1)]

	return crossedover


def recomb(surv):
	"""
	This function creates offspring through pairing of parents (diploid) and recombination (i.e, meiosis)
	"""
	# crossing over within diploid parents to make gametes, input is 1000 rows, output also 1000 rows; n_loci columns
	chroms_stacked = np.concatenate(np.stack((surv[0], surv[1]), axis = 1)).reshape((np.shape(surv[0])[0]), 2, (np.shape(surv[0])[1])) # this places the related rows together, 1st row of each together, then 2nd, then 3rd... - 4 loci. stacks the two arrays vertically 

	#recombination without looping 
	surv_stacked1 = shuffle_along_axis(chroms_stacked, 1)
	surv_stacked = np.reshape(surv_stacked1, (np.shape(chroms_stacked)[0] * 2, np.shape(chroms_stacked)[2]))
	
	# #recombination through looping
	# c = np.random.randint(1, size = np.size(surv_stacked, 1)).reshape(1, np.size(surv_stacked, 1)) # create a random array to save the results of each loop

	# x = 0
	# while x < (N_adapts[0] * 2 + 1): 
		
	# 	b = shuffle_along_axis(surv_stacked[x:(x+2)], axis = 0) #shuffle along columns in groups of 2. each group of 2 represents chrom1 an chrom2 of each individual 
	# 	c = np.concatenate((c, b), axis = 0) #update what surv_stacked is after each loop of recombination 
	# 	x+=2

	# surv_stacked = c[1:(N_adapts[0] * 2 + 1)] #remove the empty array from the top surv_stacked, update surv_stacked accordingly. chrom1, chrom2, chrom1, chrom2 seklinde devam ediyor rowlar. 

	surv_chrom1 = surv_stacked[::2] #this selects every odd row - chrom1 of N_adapts number of individuals after shuffling, both parent1 and parent2. number of rows = N_adapts, number of columns = number of loci   
	surv_chrom2 = surv_stacked[1::2] #this selects every even row - chrom 2 of N_adapts number of individuals, both parent1 and parent2 

	surv_stacked = np.hstack((surv_chrom1, surv_chrom2)) #this horizontally places chrom1 and chrom2. each row is chrom1 and chrom2 of an individual. left part = chrom1, right part = chrom2
	
	# pick random pairs of parents 
	pairs = np.resize(np.random.choice(len(surv_stacked), size=len(surv_stacked), replace=False), (int(len(surv_stacked)/2), 2))

	#determine the gamates of parents 
	parent1_gamete1 = surv_chrom1[pairs[:, 0]] #pick the chrom1 of the parents (pick the related rows, 0th element of pairs )
	parent1_gamete2 = surv_chrom2[pairs[:, 0]]

	parent2_gamete1 = surv_chrom1[pairs[:, 1]]
	parent2_gamete2 = surv_chrom2[pairs[:, 1]]

	parent1 = np.vstack((parent1_gamete1, parent1_gamete2)) #vertically add the two gametes. this gives overall genotype of parents. each row is one gamete 
	parent2 = np.vstack((parent2_gamete1, parent2_gamete2))

	#from which parent offsprings inherit each allele
	off1_chroms = np.random.randint(2, size=(len(pairs), 2)) # from which chrom to draw for off1 [p1, p2] #gives an array of 1 and 0. number of rows = number of rows of pairs. number of columns = number of rows of surv_stacked
	off2_chroms = abs(1-off1_chroms) #opposite of rand for the other offspring (the other offspring inherits the alleles from the other parent for the related loci)

	off_chrom1 = np.stack((off1_chroms[:, 0], off2_chroms[:, 0]), axis = 1).reshape(N_adapts[0], 1)
	off_chrom2 = np.stack((off1_chroms[:, 1], off2_chroms[:, 1]), axis = 1).reshape(N_adapts[0], 1)

	#create the related indices 
	even_nums = np.repeat(np.arange(0, (N_adapts[0] - 1), 2), 2).reshape(N_adapts[0], 1) #produce a list of even numbers from 0 to Nadapts - 1, not including the stop. 

	off_chrom1_index = off_chrom1 + even_nums.reshape(N_adapts[0], 1)
	off_chrom2_index = off_chrom2 + even_nums.reshape(N_adapts[0], 1)

	off = np.hstack((parent1[(off_chrom1_index)], parent2[(off_chrom2_index)])).reshape(N_adapts[0], np.shape(parent1_gamete1)[1] * 2)       #(np.size(off1_chroms, 1))) #stack the same rows from two arrays together and reformat. each row is one offspring. 
	
	return off


	# ORIGINAL CODE: 
	# pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	# rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	# rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	# off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	# off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	# off = np.append(off_1, off_2, axis=0) #each product of meiosis, diploid offspring
	# return off


def mutate(off, u_adapt, alpha, n, mut):
	"""
	This function creates mutations and updates population
	"""
	rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring [or loci??]. (creates number of off random numbers as between 0 and 1) size = number of columns of off (number of individuals)
	nmuts = np.sum(rand3 < u_adapt) # mutate if random number is below mutation rate; returns number of new mutations. 
	whomuts = np.where(rand3 < u_adapt) #indices of mutants. each refer to the index of the individuals in the off matrix. 
	newmuts = np.random.normal(0, alpha, size = (nmuts, n)) #phenotypic effect of new mutations. 0=mean, alpha=sd (how far you from the mean) rows:nmuts coloumns=n. each pair is x and y coordinates. they tell you how far and which direction you go away from the origin 
	
	#pop = np.append(off, np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants. from the identity array, only pick the rows of individuls that had a lower nmuts value than alpha. then append them next to the individuals themselves. 
	pop_chrom1 = np.split(off, 2, axis = 1)[0] #split the off array into 2 column wise. left side chrom1, right is chrom2. 
	pop_chrom2 = np.split(off, 2, axis = 1)[1]
	
	#update pop_chrom1
	added_muts = np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]) #pick the rows of mutated individuals from the identity array
	pop_chrom1 = np.append(pop_chrom1, added_muts, axis=1) #update pop_chrom1 by appending the mutation matrix

	#update pop_chrom2
	zero = np.zeros(N_adapts[0] * np.shape(added_muts)[1]).reshape(N_adapts[0], np.shape(added_muts)[1]).astype(int) #create an array of zeros. rows: n_adapts columns: same as added_muts. chrom2 doesnt mutate, so we add the zeros array. 
	pop_chrom2 = np.append(pop_chrom2, zero, axis = 1) #append zero array to chrom2

	#append pop_chrom1 and pop_chrom2 horizontally to make the pop matrix. each row is one individual. each row has the both chromosomes. left: chrom1 right: chrom2
	pop_genotype = np.append(pop_chrom1, pop_chrom2, axis = 1) #both chromosomes present

	pop_overall = (pop_chrom1 + pop_chrom2) / 2 #chromosomes averaged 

	mut = np.vstack((mut, newmuts)) #append effect of new mutations to mutation list
	
	return [newmuts, pop_genotype, pop_overall, mut]

def which_index(pop):
	return np.array([
		i for i in range(len(pop))
		if pop[i] == False ])

def add_h(pop_overall, pop_h, h):
 
	if h == "variable" or "options":
		# make a new array where all 0.5 are 0. 
		pop_overall_zero = np.copy(pop_overall)
		pop_overall_zero[pop_overall_zero == 0.5] = 0 

		# make a new array where all 1 are 0. we will sum these later to have the overall again. 
		pop_overall_h = np.copy(pop_overall) 
		pop_overall_h[pop_overall_h == 1] = 0 

		for x in range(0, np.shape(pop_overall_h)[1]):
			pop_overall_h[:, x - 1][pop_overall_h[:, x - 1] == 0.5] = pop_h[x - 1]

		# pop_overall_summed is the same as pop_overall but only h values added 
		pop_overall_summed = pop_overall_h + pop_overall_zero

	else:
		pop_overall_summed = np.copy(pop_overall)
		pop_overall_summed[pop_overall_summed == 0.5] = h

	return [pop_overall_summed]


def remove_muts(pop_genotype, mut): #here pop is the same thing as off
	"""
	This function creates mutations and updates population
	"""
	#convert mut from a list to an array 
	mut = np.array(mut)

	#split the pop into chrom1 and chrom2. 
	pop_chrom1 = np.split(pop_genotype, 2, axis = 1)[0]
	pop_chrom2 = np.split(pop_genotype, 2, axis = 1)[1]

	chrom1_zero = np.where(~pop_chrom1.any(axis = 0))[0] #find the columns that are all zeros
	chrom2_zero = np.where(~pop_chrom2.any(axis = 0))[0]

	remove = np.intersect1d(chrom1_zero, chrom2_zero) #find the indices where the two arrays have 0 in the same columns 

	if len(remove) == 0:
		pop_genotype = [pop_chrom1, pop_chrom2]

	else: 
		pop_chrom1 = np.delete(pop_chrom1, remove, 1) # delete the columns that are all zeros 
		pop_chrom2 = np.delete(pop_chrom2, remove, 1)
		pop_genotype = [pop_chrom1, pop_chrom2] #horizontally reattach the chromosomes and make the pop

		mut = np.delete(mut, remove, 0) #remove the lost loci by removing the related rows

	pop_overall = (pop_chrom1 + pop_chrom2) / 2 
	
	return[pop_chrom1, pop_chrom2, pop_genotype, pop_overall, mut, remove]

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

nreps = 1 #number of replicates for each set of parameters
ns = [2] #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

N_adapts = [1000] #number of diploid individuals (positive integer)
alpha_adapts = [0.10] #mutational sd (positive real number)
u_adapts = [0.001] #mutation probability per generation per genome (0<u<1). if this is 0.5, this means half of the population is likely to mutate 
sigma_adapts = [1, 3, 5] #selection strengths

opt_dists = [1] #distance to optima

n_angles = 12 #number of angles between optima to simulate (including 0 and 180) (>=2)

maxgen = 1400 #total number of generations populations adapt for
 
# 9 -> variable, chosen from random distribution
# options -> 0, 0.5, 1
#
dom = [9 , 0.5]

######################################################################
##PARAMETERS FOR HYBRIDS##
######################################################################

nHybrids = 200 #number of hybrids to make at end of each replicate

######################################################################
##FUNCTION FOR POPULATIONS TO ADAPT##
######################################################################

def main():

	#open the files, return file handles to each 
	[fileHandle_A, fileHandle_B, fileHandle_C, fileHandle_D, fileHandle_E] = open_output_files(ns, N_adapts, alpha_adapts, u_adapts, sigma_adapts, opt_dists, dom)

	#loop over population size
	i_N = 0
	while i_N < len(N_adapts):
		N_adapt = N_adapts[i_N]

		#loop over h values
		i_dom = 0 
		while i_dom < len(dom):
			h = dom[i_dom]

			#loop over selection strength
			i_sigma = 0
			while i_sigma < len(sigma_adapts):
				sigma_adapt = sigma_adapts[i_sigma]

				#loop over mutation rate (u)
				i_u = 0 
				while i_u < len(u_adapts):
					u_adapt = u_adapts[i_u]

					#loop over mutational sd
					i_alpha = 0 
					while i_alpha < len(alpha_adapts):
						alpha_adapt = alpha_adapts[i_alpha]

						#loop over opt_dist
						i_opt_dist = 0 
						while i_opt_dist < len(opt_dists):
							opt_dist = opt_dists[i_opt_dist]

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

										mut = np.array([[0] * n]) #similar to above. filled with zeroes. number of coloumns: n. rows: 1. convert to a list 

										[pop1, mut1] = [popfound, mut] # this creates pop and mut arrays for both parents. they are the same because we start from the same point. 
										[pop2, mut2] = [popfound, mut] # mut1 = how farther you go from the origin due to mutations in pop1. same for mut2

										pop1_chrom1 = popfound1 # genotype of 1st chromosome of pop1
										pop1_chrom2 = popfound1 # genotype of 2nd chromosome of pop1

										pop2_chrom1 = popfound2 # 1st chromosome of pop2
										pop2_chrom2 = popfound2 # 2nd chromosome of pop2

										pop1_overall = ((pop1_chrom1 + pop1_chrom2) / 2 ) # two chromosomes of pop1 averaged
										pop2_overall = ((pop2_chrom1 + pop2_chrom2) / 2 ) # two chromosomes of pop2 averaged

										pop1_overall_summed = np.copy(pop1_overall)
										pop2_overall_summed = np.copy(pop2_overall)

										if h == 9:
											# this is the only time we use these first values 
											pop1_first_h = np.round(np.random.uniform(low = 0, high = 1, size = 1), 2)
											pop2_first_h = np.round(np.random.uniform(low = 0, high = 1, size = 1), 2)

											pop1_h = np.copy(pop1_first_h)
											pop2_h = np.copy(pop2_first_h)

										elif h == 'options': 
											options = np.array([0, 1, 0.5])

											pop1_first_h = np.round(np.random.choice(options, 1), 2)
											pop2_first_h = np.round(np.random.choice(options, 1), 2)

											pop1_h = np.copy(pop1_first_h)
											pop2_h = np.copy(pop2_first_h)

										else: 
											# h is an integer and not picked from random distribution. attach the first h value to the pop_h matrices 
											pop1_h = h
											pop2_h = h 

										# intitialize generation counter
										gen = 0

										# run until maxgen
										while gen < maxgen + 1:

											# genotype to phenotype (diploid):
											phenos1 = np.dot(pop1_overall_summed, mut1) #sum mutations held by each individual
											phenos2 = np.dot(pop2_overall_summed, mut2) #sum mutations held by each individual

											# phenotype to fitness
											w1 = fitness(phenos1, theta1, sigma_adapt)
											w2 = fitness(phenos2, theta2, sigma_adapt)

											# wright-fisher (multinomial) sampling
											# number of times each parent chosen, drawing samples from a multinomial ditribution
											# N_adapt = number of experiments, w1/sum(w1 = probability of parent1 being chosen. if you are more fit, you are chosen more often. 
											parents1 = np.random.multinomial(N_adapt, w1/sum(w1))
											off1_chrom1 = np.repeat(pop1_chrom1, parents1, axis = 0) 
											off1_chrom2 = np.repeat(pop1_chrom2, parents1, axis = 0)
											off1 = [off1_chrom1, off1_chrom2]
											
											parents2 = np.random.multinomial(N_adapt, w2/sum(w2)) # number of times each parent chosen
											off2_chrom1 = np.repeat(pop2_chrom1, parents2, axis = 0) # offspring genotypes of pop2
											off2_chrom2 = np.repeat(pop2_chrom2, parents2, axis = 0)
											off2 = [off2_chrom1, off2_chrom2]

											# mating and recombination
											off1 = recomb(off1)
											off2 = recomb(off2)

											# mutation and population update
											[pop1_newmuts, pop1_genotype, pop1_overall, mut1] = mutate(off1, u_adapt, alpha_adapt, n, mut1)
											[pop2_newmuts, pop2_genotype, pop2_overall, mut2] = mutate(off2, u_adapt, alpha_adapt, n, mut2)

											if h == 9:
												# generate new h values
												pop1_h_value = np.round(np.random.uniform(low = 0, high = 1, size = len(pop1_newmuts + 1)).reshape(len(pop1_newmuts), 1), 2)
												pop2_h_value = np.round(np.random.uniform(low = 0, high = 1, size = len(pop2_newmuts + 1)).reshape(len(pop2_newmuts), 1), 2)

												# this the overall list where we save all the h values 
												pop1_h = np.append(pop1_h, pop1_h_value)
												pop1_h = np.reshape(pop1_h, (np.shape(pop1_h)[0], 1)) #reshape into one column 
												
												pop2_h = np.append(pop2_h, pop2_h_value)
												pop2_h = np.reshape(pop2_h, (np.shape(pop2_h)[0], 1)) #reshape into one column 

											elif h == 'options': 
												# generate new h values
												pop1_h_value = np.round(np.random.choice(options, len(pop1_newmuts)).reshape(len(pop1_newmuts), 1), 2)
												pop2_h_value = np.round(np.random.choice(options, len(pop2_newmuts)).reshape(len(pop2_newmuts), 1), 2)

												# this the overall list where we save all the h values 
												pop1_h = np.append(pop1_h, pop1_h_value)
												pop1_h = np.reshape(pop1_h, (np.shape(pop1_h)[0], 1)) #reshape into one column 
												
												pop2_h = np.append(pop2_h, pop2_h_value)
												pop2_h = np.reshape(pop2_h, (np.shape(pop2_h)[0], 1)) #reshape into one column 


											else:
												pop1_h = np.append(pop1_h, np.repeat(h, len(pop1_newmuts))) 
												pop1_h = np.reshape(pop1_h, (np.shape(pop1_h)[0], 1)) #reshape into one column

												pop2_h = np.append(pop2_h, np.repeat(h, len(pop2_newmuts)))
												pop2_h = np.reshape(pop2_h, (np.shape(pop2_h)[0], 1)) #reshape into one column
											
											# remove lost mutations (all zero columns in pop)
											[pop1_chrom1, pop1_chrom2, pop1_genotype, pop1_overall, mut1, remove1] = remove_muts(pop1_genotype, mut1)
											[pop2_chrom1, pop2_chrom2, pop2_genotype, pop2_overall, mut2, remove2] = remove_muts(pop2_genotype, mut2)
											
											[pop1_overall_summed] = add_h(pop1_overall, pop1_h, h)
											[pop2_overall_summed] = add_h(pop2_overall, pop2_h, h)

											pop1_h = np.delete(pop1_h, remove1, 0) #remove the same rows from the pop_h matrix 
											pop2_h = np.delete(pop2_h, remove2, 0)

											# go to next generation
											gen += 1

											#############################################
											# CREATE THE ADAPTATION SUMMARY (of parents)
											#############################################
											
											one = np.ones(len(mut1)).reshape(len(mut1), 1)
											two = np.random.randint(low = 2, high = 3, size = len(mut2)).reshape(len(mut2), 1)
											pop1_or_pop2 = np.vstack((one, two))

											pop1_freq = (np.mean(pop1_overall, axis = 0)).reshape(np.shape(pop1_overall)[1], 1)
											pop2_freq = (np.mean(pop2_overall, axis = 0)).reshape(np.shape(pop2_overall)[1], 1)
											freq_overall = np.round(np.vstack((pop1_freq, pop2_freq)), 4)
											
											pop1_mut_n1 = mut1[:, 0].reshape(np.shape(mut1)[0], 1)
											pop2_mut_n1 = mut2[:, 0].reshape(np.shape(mut2)[0], 1)
											mut_n1_overall = np.round(np.vstack((pop1_mut_n1, pop2_mut_n1)), 3)
											
											pop1_mut_n2 = mut1[:, 1].reshape(np.shape(mut1)[0], 1)
											pop2_mut_n2 = mut2[:, 1].reshape(np.shape(mut2)[0], 1)
											mut_n2_overall = np.round(np.vstack((pop1_mut_n2, pop2_mut_n2)), 3)

											h_overall = np.vstack((pop1_h, pop2_h))

											freq_fixed = np.where(freq_overall > 0.5)[0]
											h_fixed = h_overall[freq_fixed]
											h_fixed_mean = np.mean(h_fixed)
											
											
											if h == 0:
												adapt_summary_0 = np.column_stack((pop1_or_pop2, freq_overall, mut_n1_overall, mut_n2_overall, h_overall))

											elif h == 1:
												adapt_summary_1 = np.column_stack((pop1_or_pop2, freq_overall, mut_n1_overall, mut_n2_overall, h_overall))

											elif h == 0.5:
												adapt_summary_half = np.column_stack((pop1_or_pop2, freq_overall, mut_n1_overall, mut_n2_overall, h_overall))

											else:
												adapt_summary_variable = np.column_stack((pop1_or_pop2, freq_overall, mut_n1_overall, mut_n2_overall, h_overall))

										#HYBRIDIZATION - always choose parents from different populations 

										#genotype of hybrids
										pop1_chrom1 = pop1_genotype[0]
										pop1_chrom2 = pop1_genotype[1]

										pop2_chrom1 = pop2_genotype[0]
										pop2_chrom2 = pop2_genotype[1]
										
										#make the zero matrices
										pop1_zero1 = np.zeros(len(pop1_chrom1) * pop2_chrom1.shape[1]).reshape(len(pop2_chrom1), pop2_chrom1.shape[1])
										pop1_zero2 = np.zeros(len(pop1_chrom2) * pop2_chrom2.shape[1]).reshape(len(pop2_chrom2), pop2_chrom2.shape[1])

										pop2_zero1 = np.zeros(len(pop2_chrom1) * pop1_chrom1.shape[1]).reshape(len(pop1_chrom1), pop1_chrom1.shape[1])
										pop2_zero2 = np.zeros(len(pop2_chrom2) * pop1_chrom2.shape[1]).reshape(len(pop1_chrom2), pop1_chrom2.shape[1])


										#attach the zero matrices
										pop1_chrom1_has0 = np.hstack((pop1_chrom1, pop1_zero1))
										pop1_chrom2_has0 = np.hstack((pop1_chrom2, pop1_zero2))

										pop2_chrom1_has0 = np.hstack((pop2_zero1, pop2_chrom1))
										pop2_chrom2_has0 = np.hstack((pop2_zero2, pop2_chrom2))

										#make pairs
										pairs_hybrid = np.resize(np.random.choice(len(pop1_chrom1), size=len(pop1_chrom1), replace=False), (int(len(pop1_chrom1)/2), 2))

										#pick the related chromosomes of the pairs  
										pop1_chrom1_hy = pop1_chrom1_has0[pairs_hybrid[:, 0]]
										pop1_chrom2_hy = pop1_chrom2_has0[pairs_hybrid[:, 0]]

										pop2_chrom1_hy = pop2_chrom1_has0[pairs_hybrid[:, 1]]
										pop2_chrom2_hy = pop2_chrom2_has0[pairs_hybrid[:, 1]]

										#recombination:
										#randomly pick 0 or 1 to decide which pairs to match 

										num = np.random.randint(2, size = 1).tolist()

										if num[0] == 0:
											F1_before_recomb1 = np.concatenate(np.stack((pop1_chrom1_hy, pop2_chrom1_hy), axis = 1))
											F1_after_recomb1 = crossover(F1_before_recomb1)
											F1_after_recomb1_chrom1 = F1_after_recomb1[::2] #picks every other odd row, chrom1
											F1_after_recomb1_chrom2 = F1_after_recomb1[1::2] #picks every other even row, chrom2

											
											F1_before_recomb2 = np.concatenate(np.stack((pop1_chrom2_hy, pop2_chrom2_hy), axis = 1))
											F1_after_recomb2 = crossover(F1_before_recomb2)
											F1_after_recomb2_chrom1 = F1_after_recomb2[::2] #picks every other odd row, chrom1
											F1_after_recomb2_chrom2 = F1_after_recomb2[1::2] #picks every other even row, chrom2


										elif num[0] == 1:
											F1_before_recomb1 = np.concatenate(np.stack((pop1_chrom1_hy, pop2_chrom2_hy), axis = 1))
											F1_after_recomb1 = crossover(F1_before_recomb1)
											F1_after_recomb1_chrom1 = F1_after_recomb1[::2] #picks every other odd row, chrom1
											F1_after_recomb1_chrom2 = F1_after_recomb1[1::2] #picks every other even row, chrom2

											F1_before_recomb2 = np.concatenate(np.stack((pop1_chrom2_hy, pop2_chrom1_hy), axis = 1))
											F1_after_recomb2 = crossover(F1_before_recomb2)
											F1_after_recomb2_chrom1 = F1_after_recomb2[::2] #picks every other odd row, chrom1
											F1_after_recomb2_chrom2 = F1_after_recomb2[1::2] #picks every other even row, chrom2

										#hybrid phenotypes
										mut_hybrid = np.vstack((mut1, mut2))

										#averaged hybrid genotypes to calculate the hybrid phenotypes
										hybrid_overall_recomb1 = (F1_after_recomb1_chrom2 + F1_after_recomb1_chrom1) / 2
										hybrid_overall_recomb2 = (F1_after_recomb2_chrom2 + F1_after_recomb2_chrom1) / 2

										hybrid_overall_all = np.vstack((hybrid_overall_recomb1, hybrid_overall_recomb2)) #stack all hybrids phenos. each row is one individual. 

										# generate the hybrid h values matrix. 
										hybrid_h = np.vstack((pop1_h, pop2_h))

										# split the overall into 2 chromosomes accordingly
										hybrid_overall_one = hybrid_overall_all[:, 0:(np.shape(pop1_h)[0])]
										hybrid_overall_two = hybrid_overall_all[:, (np.shape(pop1_h)[0]):((np.shape(pop1_h)[0]) + (np.shape(pop2_h)[0]))]

										# replace with the h values 
										for x in range(0, np.shape(hybrid_overall_one)[1]):
											hybrid_overall_one[:, x - 1][hybrid_overall_one[:, x - 1] == 0.5] = pop1_h[x - 1]

										for x in range(0, np.shape(hybrid_overall_two)[1]):
											hybrid_overall_two[:, x - 1][hybrid_overall_two[:, x - 1] == 0.5] = pop2_h[x - 1]

										# connect the two chromosomes again to make the overall 
										hybrid_overall_all = np.hstack((hybrid_overall_one, hybrid_overall_two))

										hybrid_pheno = np.dot(hybrid_overall_all, mut_hybrid)

										# find the mean phenos of each column 
										phenos1_1 = np.array([np.round((np.sum(phenos1, axis = 0)[0]) / len(phenos1), 3)])
										phenos1_2 = np.array([np.round((np.sum(phenos1, axis = 0)[1]) / len(phenos1), 3)])

										phenos2_1 = np.array([np.round((np.sum(phenos2, axis = 0)[0]) / len(phenos2), 3)])
										phenos2_2 = np.array([np.round((np.sum(phenos2, axis = 0)[1]) / len(phenos2), 3)])

										hybrid_phenos1 = np.array([np.round((np.sum(hybrid_pheno, axis = 0)[0]) / len(hybrid_pheno), 3)])
										hybrid_phenos2 = np.array([np.round((np.sum(hybrid_pheno, axis = 0)[1]) / len(hybrid_pheno), 3)])

										# when h is variable, find the average h values 
										h_averaged = (np.mean(h_overall, axis = 0))

										# HYBRID FITNESS 
										p_fit1 = w1
										p_fit2 = w2 

										# calculate the hybrid fitness in two optima and then assign them the higher of the two values 
										h_fit1 = np.array(fitness(hybrid_pheno, theta1, sigma_adapt))
										h_fit2 = np.array(fitness(hybrid_pheno, theta2, sigma_adapt))
										
										# make a random array as long as the h_fit1. We will replace the numbers here as we run the loop below. 
										h_fit = np.ones(len(h_fit1))

										# go through the values in both h_fit1 and h_fit2. Pick the greater value. Save the values in the h_fit array. 
										i = 0 
										while i < len(h_fit1):
											if h_fit1[i] > h_fit2[i]:
												h_fit[i] = h_fit1[i]
											else: 
												h_fit[i] = h_fit2[i]
											i += 1

										# reshape the fitness arrays. only 1 column. 
										h_fit = np.reshape(h_fit, (len(h_fit), 1))
										w1 = np.reshape(w1, (len(w1), 1))
										w2 = np.reshape(w2, (len(w2), 1))

										# SAVE PHENO DATA / PHENOTYPES 
										pheno_data = np.column_stack([np.array([i+1 for i in range(len(hybrid_pheno))]), np.array([rep+1]*len(hybrid_pheno)), np.array([round(angles[j]*180/math.pi,2)]*len(hybrid_pheno)), np.array([opt_dist * (2*(1-math.cos(angles[j])))**(0.5)]*len(hybrid_pheno)), hybrid_pheno])
										np.savetxt(fileHandle_B, pheno_data, fmt = '%.3f', delimiter = ',')
										
										# SAVE THE H VALUES / SUMMARY 
										sum_data = np.column_stack(np.array(([round(angles[j]*180/math.pi,2), rep+1, np.round(opt_dist * (2*(1-math.cos(angles[j])))**(0.5), 3), phenos1_1, phenos1_2, phenos2_1, phenos2_2, hybrid_phenos1, hybrid_phenos2, np.round(h_averaged, 3), h_fixed_mean, u_adapt, sigma_adapt, alpha_adapt, opt_dist])))
										np.savetxt(fileHandle_A, sum_data, fmt = '%.3f', delimiter = ',') 
										
										# SAVE THE MUT SUMMARY / ADAPTATION
										adapt_summary = np.vstack((adapt_summary_variable))
										np.savetxt(fileHandle_C, adapt_summary, fmt = '%.3f', delimiter = ',')

										#SAVE THE FITNESS VALUES 
										#put the other data in the correct format
										a = np.reshape(np.repeat(N_adapt, np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)
										b = np.reshape(np.repeat(sigma_adapt, np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)
										c = np.reshape(np.repeat(u_adapt, np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)
										d = np.reshape(np.repeat(alpha_adapt, np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)
										e = np.reshape(np.repeat(opt_dist, np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)
										f = np.reshape(np.repeat(round(angles[j]*180/math.pi,2), np.shape(w1)[0], axis = 0), np.shape(w1)[0], 1)

										fitness_sum = np.column_stack((w1, w2, h_fit, f, a, b, c, d, e))
										np.savetxt(fileHandle_D, fitness_sum, fmt = '%.3f', delimiter = ',')

										#CREATE THE FIT SUM ARRAY 
										#calculate the mean fitnesses 
										w1_mean = np.mean(fitness_sum[:, 0].reshape(np.shape(fitness_sum)[0], 1), axis = 0)
										w2_mean = np.mean(fitness_sum[:, 1].reshape(np.shape(fitness_sum)[0], 1), axis = 0)
										h_fit_mean = np.mean(fitness_sum[:, 2].reshape(np.shape(fitness_sum)[0], 1), axis = 0)

										fit_mean = np.column_stack((w1_mean, w2_mean, h_fit_mean, h, np.round(angles[j]*180/math.pi,2), N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n)) # , (round(angles[j]*180/math.pi,2), N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n)))
										np.savetxt(fileHandle_E, fit_mean, fmt = '%.3f', delimiter = ',')

										# print an update
										if h == 9:
											print('N = %d, sigma = %.3f, u = %.3f, alpha = %.3f, opt_dist = %.2f, n=%d, angle=%r, rep=%d, h= %d' %(N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n, round(angles[j]*180/math.pi,2), rep+1, h)) 
										else:
											print('N = %d, sigma = %.3f, u = %.3f, alpha = %.3f, opt_dist = %.2f, n=%d, angle=%r, rep=%d, h= %s' %(N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n, round(angles[j]*180/math.pi,2), rep+1, h)) 

										# go to next rep
										rep += 1

									#go to next optimum
									j += 1

								#next dimension
								l += 1

							#next opt_dist
							i_opt_dist += 1 
						
						#next alpha_adapt value 
						i_alpha += 1	
						
					#next u_adapt value 
					i_u += 1 

				#go to next sigma value
				i_sigma += 1

			#go to the next h value 
			i_dom += 1

		#go to next N value
		i_N += 1

		#clean up
		close_output_files(fileHandle_A)
		close_output_files(fileHandle_B) 
		close_output_files(fileHandle_C)
		close_output_files(fileHandle_D)
		close_output_files(fileHandle_E)

		# add header to the files (column names)
		os.chdir('/Users/student/Desktop/asli/summer-2019/data_done')

		#summary
		sim_id = '_n%s_N%s_alpha%s_u%s_sigma%s_opt_dist%s_dom%s' %(ns, N_adapts, alpha_adapts, u_adapts, sigma_adapts, opt_dists, dom)
		df_summary = pd.read_csv("mutsummary%s.csv" %(sim_id), header = None)
		df_summary.columns = ['angle', 'reps', 'opt_dist', 'phenos1_1', 'phenos1_2', 'phenos2_1', 'phenos2_2', 'hybrid_phenos1', 'hybrid_phenos2', 'h_mean', 'h_fixed_mean', 'u', 'sigma', 'alpha', 'opt_dist'] 
		df_summary.to_csv("mutsummary%s.csv" %(sim_id))

		#adaptation
		df_adapt = pd.read_csv("adaptation%s.csv" %(sim_id), header = None)
		df_adapt.columns = ['pop1_or_pop2', 'freq_overall', 'mut_n1_overall', 'mut_n2_overall', 'h_overall'] 
		df_adapt.to_csv("adaptation%s.csv" %(sim_id))

		#phenotypes
		df_pheno = pd.read_csv("phenotypes%s.csv" %(sim_id), header = None)
		df_pheno.columns = ['individuals', 'reps', 'angle', 'opt_dist', 'hybrid_phenos1', 'hybrid_phenos2']
		df_pheno.to_csv("phenotypes%s.csv" %(sim_id))

		#individual fitness
		df_fitness = pd.read_csv("fitness%s.csv" %(sim_id), header = None)
		df_fitness.columns = ['parent1', 'parent2', 'hybrids', 'angle', 'N_adapt', 'sigma', 'u', 'alpha', 'opt_dist']
		df_fitness.to_csv("fitness%s.csv" %(sim_id))

		#fitMean
		df_fitMean = pd.read_csv("fitMean%s.csv" %(sim_id), header = None)
		df_fitMean.columns = ['parent1', 'parent2', 'hybrids', 'h', 'angle', 'N_adapt', 'sigma_adapt', 'u_adapt', 'alpha_adapt', 'opt_dist', 'n']
		df_fitMean.to_csv("fitMean%s.csv" %(sim_id))



		
######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    

start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))