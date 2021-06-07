import numpy as np
import time
import csv
import math
import pandas as pd
import os
from numpy.random import RandomState
import sys
from scipy.integrate import quad
import random

# REPLACE THIS: 
data_dir = '/Users/asli_munzur/Desktop/BIOL448/data'

######################################################################
## PREPARE THE DIRS ##
######################################################################

# REMOVE EMPTY SUBDIRS IN DIRECTORY ABOVE  
folders = list(os.walk(data_dir))[1:]

for folder in folders:
	# folder example: ('FOLDER/3', [], ['file'])
	if not folder[2]:
		os.rmdir(folder[0])

# go to this dir and count the number FOLDERS in there 
# this will help determine how many runs have been made before so that we dont overwrite any of the previous data 
files = folders = 0

for _, dirnames, filenames in os.walk(data_dir):
	# ^ this idiom means "we won't be using this value"
	files += len(filenames)
	folders += len(dirnames)

trial = "trial" + str(folders) # now based on number of already existing folders, figure out the trail number
updated_path = data_dir + "/" + trial

# if the path we defined above exists, increase the trial number one by one until you find a unique path number 
while os.path.exists(updated_path): # while true
	folders = folders + 1 
	trial = "trial" + str(folders)
	updated_path = data_dir + "/" + trial

	if os.path.exists(updated_path) == False: 
		break

os.makedirs(updated_path)

######################################################################
## FUNCTIONS ##
######################################################################


def open_output_files(data_dir, n, N, alpha, sigma, opt_dist, dom, pevo):
	"""
	This function opens the output files and returns file
	handles to each.
	"""

	# update here if directory changes 
	data_dir = data_dir

	outfile_A = open(os.path.join(updated_path, "mutsummary.csv"), "w") #summary data
	outfile_B = open(os.path.join(updated_path, "phenotypes.csv"), "w") #hybrid phenotypes at the end 
	outfile_C = open(os.path.join(updated_path, "adaptation.csv"), "w")
	outfile_D = open(os.path.join(updated_path, "F2_BC.csv"), "w")#save the fitness of the parents and the F1 hybrids in the same file 
	outfile_E = open(os.path.join(updated_path, "fitMean.csv"), "w")

	return [outfile_A, outfile_B, outfile_C, outfile_D, outfile_E] #later, add the outfile D if needed

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

def fitness_F2(phenos, theta, sigma):
	"""
	This function determines relative fitness
	"""
	dist = np.linalg.norm(phenos - theta) #phenotypic distance from optimum
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
	surv is the genotype, given in a list. surv[0] is the 1st kromozom, surv[1] is the 2nd kromozom. 
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

def mutate(alpha, n, dominance):  
	"""
	# generates same mutations in the same order in both populations, depending on what 
	"""
	newmut = np.random.normal(0, alpha, size = (1, n)) # phenotypic effect of new mutations, pick as many as n
	newmut.shape = (1, n) # reshape into intp 1 row n cols

	if dominance == 0.5:
		h_z = 0.5
	elif dominance == 8: 
		h_z = np.random.uniform(0.65, 1, size = 1)
	elif dominance == 9: 
		h_z = np.random.uniform(0, 1, size = 1)
		# h_z = 0.5
	
	return(h_z, newmut) # return the h value and phenotypic value for the new mutation 


def compute_s(newmut, WT_pheno, theta, sigma):
	
	# call fitness function on both phenos
	mutant_pheno = WT_pheno + newmut

	w_homo_mutant = fitness(mutant_pheno, theta, sigma)
	w_WT = fitness(WT_pheno, theta, sigma)
	
	s = (w_homo_mutant / w_WT) - 1
	
	return s

def compute_hfit(mutant_pheno, WT_pheno, theta, sigma, s, h_z):
	"""
	# this function computes the h_fit values based on fitness, s and h_z values
	"""
	
	# compute the new phenotype using wildtype and mutant phenotypes
	new_pheno = WT_pheno + h_z * mutant_pheno
	
	# call fitness function on both phenos
	w_mutant = fitness(new_pheno, theta, sigma)
	w_WT = fitness(WT_pheno, theta, sigma)
	
	h_fit = ((w_mutant / w_WT) - 1) / s
	
	return(h_fit)

# h_fit is the value that we compute, not the value we randomly generate for each mut (h_z)
def did_it_fix_OLD(N, s, h_fit):
	
	"""
	# this function decides if a mutation fixes or not, based on a couple of parameters 
	# old equation, we no longer use it
	"""
	
	# define the function we want to integrate, part of the intense function
	def f(x, N, s, h_fit):
		return math.e**(-2*N*s*(2*h_fit-1)*x*(1-x)-(2*N*s*x))
	
	p = 1/(2*N)

	# call quad to integrate the function
	res1 = quad(f, 0, p, args = (N, s, h_fit))[0] # numerator
	res2 = quad(f, 0, 1, args = (N, s, h_fit))[0] # denominator

	# divide the results of the two integrations to get the fixation probability
	result = res1/res2
 
	rand = np.random.uniform(size = 1) # pick a random number between 0, 1 - random uniform
	if rand < result: # mutation fixes here if the random number is more than the probability we computed
		did_it_fix = int(1) # if it fixes, give it 1
	  
	else: # mutation doesn't fix
		did_it_fix = int(0) # if it gets lost, give it 0
		
	return did_it_fix

def did_it_fix(h_fit, s):

	result = 2*h_fit*s

	rand = np.random.uniform(size = 1) # pick a random number between 0, 1 - random uniform
	if rand < result: # mutation fixes here if the random number is more than the probability we computed
		did_it_fix = int(1) # if it fixes, give it 1
	  
	else: 
		did_it_fix = int(0) # if it gets lost, give it 0

	return did_it_fix

def evaluate_muts(did_it_fix, mut, newmut, pop_h, h_z):
	
	"""
	# this function modifies the genotype and mut arrays if the mutation fixes, otherwise does nothing
	"""
	
	if did_it_fix == 1: 
		mut = np.vstack((mut, newmut)) # add the new mut to the array        
		pop_h = np.vstack((pop_h, h_z)) # add the h_z
		
	else: 
		pass # do nothing
	
	return [mut, pop_h]

def which_index(pop):
	return np.array([
		i for i in range(len(pop))
		if pop[i] == False ])

# this function is for adding the effect of h values into the pop arrays. 
def add_h(pop_overall, pop_h, h):
 
	if h == 9 or "options":
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

def dealwithh(pop_overall, pop_h):

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

	return pop_h, pop_overall_summed

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

def calc_kenmet(pop1_genotype_pe, pop1_pe_idx_checked, pop2_genotype_pe, pop2_pe_idx_checked):
	
	# pop1_shared = len(pop1_pe_idx) - 1
	pop1_total = np.shape(pop1_genotype_pe)[1]

	# pop2_shared = len(pop2_pe_idx) - 1
	pop2_total = np.shape(pop2_genotype_pe)[1]

	if pop1_total == 0 or pop2_total == 0: # if the pops didn't fix any muts automatically assign 0 to kenmet
		kenmet = 0
	else:
		kenmet = round((0.5 * (pop1_pe_idx_checked / pop1_total + pop2_pe_idx_checked / pop2_total)) * 100, 1)

	return kenmet

# this computes factors for a given number 
def factors(n):
	return list(set(
		factor for i in range(1, int(n**0.5) + 1) if n % i == 0
		for factor in (i, n//i)))

def add_modularity(mut, modularity, m):

	# MODULARITY ON 
	if modularity == 1: 

		n = np.shape(mut)[1] # extract dimensions from the mut array 

		# choosing factors randomly
		# factor_list = factors(n)
		# module_number = random.choice(factor_list) # randomly pick an element from the factors, this is how many modules we will have 

		# factor given 
		module_number = m

		zero_array = np.zeros_like(mut) # make an array of zeros, same shape as the mut array

		if module_number != 1 and module_number != n: # exclude two extreme cases, 1 is no modularity, and n is complete modularity, we evaluate them below

			i = 0
			while i < np.shape(mut)[0]: # while i < less than number of mutations

				x = int(np.random.randint(low = 1, high = module_number, size = 1)) # random number between 1 and module_number, including module_number, this helps choose which module mutation affects

				beginning_trait = (x - 1)*module_number # number of traits that come before our chosen trait
				ending_trait = int(beginning_trait + (n / module_number)) # add the size of a module, convert to an integer
				
				zero_array[i, beginning_trait:ending_trait] = mut[i, beginning_trait:ending_trait] # extract the cols from the ith row and add to 0 array 

				i += 1 # go to the next row 

		elif module_number == n: # complete modularity, each trait is a module in its own 

			i = 0
			while i < np.shape(mut)[0]: # while i < less than number of mutations

				x = int(np.random.randint(low = 0, high = module_number, size = 1)) # random number between 1 and module_number, including module_number, this helps choose which module mutation affects

				zero_array[i, x] = mut[i, x] # add the traits we extracted above to the zero array 

				i += 1 # go to the next row 


		else: # no modularity, 
			zero_array = mut

	# MODULARITY OFF 
	else: 
		zero_array = mut # do nothing, mut array stays unchanged 

	return(zero_array) # this is the updated mut array


def make_backcrosses(backcross_mut, other_mut, mut1_del, mut2_del, mut1_pe, pevo, which_parent, how_many_backcrosses, hybrid_h, mut_hybrid, theta1, theta2, sigma_adapt):
	
	# an empty list to save the genotypes
	backcross_list = list()

	i = 0 
	while i < how_many_backcrosses:
		# without pevo
		if pevo == 0: 
			parent_loci = np.shape(backcross_mut)[0]
			non_parent_loci = np.shape(other_mut)[0]

			a = np.random.choice([1, "h"], size = parent_loci, p = [0.5, 0.5])
			b = np.random.choice([0, "h"], size = non_parent_loci, p = [0.5, 0.5])

			backcross_hybrid = np.concatenate([a, b])

		# with pevo
		else: 

			a = np.ones(mut1_pe.shape[0]) # number of parallel loci

			if which_parent == "parent1":
				b = np.random.choice([1, "h"], size = mut1_del.shape[0], p = [0.5, 0.5])
				c = np.random.choice([0, "h"], size = mut2_del.shape[0], p = [0.5, 0.5])

			else:
				b = np.random.choice([1, "h"], size = mut2_del.shape[0], p = [0.5, 0.5])
				c = np.random.choice([0, "h"], size = mut1_del.shape[0], p = [0.5, 0.5])

			backcross_hybrid = np.concatenate([a, b, c])

		# add the backcross to the list 
		backcross_list.append(backcross_hybrid)
		i += 1

	# REPLACE THE "h" WITH THE ACTAUL H VALUES
	hybrid_h = hybrid_h.flatten()
	for backcross in backcross_list: 

		i = 0 
		while i < len(hybrid_h): 
			idx = np.where(backcross == "h")[0] # find the idx of "h"
			backcross[idx] = hybrid_h[idx] # replace those with the real h values

			i += 1

	backcross_list = [[float(x) for x  in backcross] for backcross in backcross_list]

	# CALCULATE PHENOTYPE
	pheno_list = []
	 
	for backcross in backcross_list:

		pheno = list(np.dot(backcross, mut_hybrid))
		pheno_list.append(pheno)

	# CALCULATE FITNESS
	fitness_list = []

	for pheno in pheno_list:

		w1 = np.array(fitness_F2(pheno, theta1, sigma_adapt))
		w2 = np.array(fitness_F2(pheno, theta2, sigma_adapt))

		# pick the higher one
		w = max(w1, w2)

		fitness_list.append(w)

	return(pheno_list, fitness_list)


######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

nreps = 30 #number of replicates for each set of parameters
ns = [2] #phenotypic dimensions (positive integer >=1)

N_adapts = [1000] #number of diploid individuals (positive integer)

alpha_adapts = [0.1] #mutational sd (positive real number)
# alpha_adapts = [0.1] #mutational sd (positive real number)

# u_adapt = (0.0001/alpha_adapt)

sigma_adapts = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10] #selection strengths

opt_dists = [1] #distance to optima

n_angles = 48 #number of angles between optima to simulate (including 0 and 180) (>=2)

dom = [0.5, 9]
# 9 -> variable, chosen from random distribution
# options -> 0, 0.5, 1

pevo = [0, 1] #this is parallel evolution. either on or off. 0 means off, 1 means on. 

modularity_all = [0]
m_all = factors(ns[0]) # all factors of n 

circling_percentage = 10
how_many_F2 = how_many_backcrosses = 500 # they must be equal

######################################################################
## MAKE THE LOG FOR THE RUN ##
######################################################################   
os.chdir(updated_path) 
f = open("log","w") # open a new text file

param_names_list = [

'nreps = ' + str(nreps), 
'ns = ' + str(ns), 
'N_adapts = ' + str(N_adapts), 
'alpha_adapts = ' + str(alpha_adapts), 
'sigma_adapts = ' + str(sigma_adapts), 
'opt_dists = ' + str(opt_dists), 
'n_angles = '+ str(n_angles), 
'dom = ' + str(dom),
'pevo = ' + str(pevo),  
'modularity = ' + str(modularity_all),
'how_many_F2 = ' + str(how_many_F2),
'file_name = ' + str(os.path.basename(sys.argv[0]))
]

for element in param_names_list:
	f.write(element + '\n')

f.close() # close the file, good practice! 

######################################################################
##FUNCTION FOR POPULATIONS TO ADAPT##
######################################################################

def main():

	#open the files, return file handles to each 
	[fileHandle_A, fileHandle_B, fileHandle_C, fileHandle_D, fileHandle_E] = open_output_files(data_dir, ns, N_adapts, alpha_adapts, sigma_adapts, opt_dists, dom, pevo)

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

				#loop over modularity
				i_modul = 0
				while i_modul < len(modularity_all):
					modularity = modularity_all[i_modul]

					#loop over mutational sd
					i_alpha = 0 
					while i_alpha < len(alpha_adapts):
						alpha_adapt = alpha_adapts[i_alpha]

						#find over mutation rate 
						u_adapt = (0.00001/alpha_adapt)
						#can uncomment this part later and add a u+=1 at the end to loop over u values 
						# i_u = 0 
						# while i_u < len(u_adapts):
						# 	u_adapt = u_adapts[i_u]

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

									#loop over all replicates
									rep = 0

									while rep < nreps:

										#loop over modularity factors
										i_m = 0
										if modularity == 0: # if there is no modularity, don't loop over the factors 
											m_all = [1]

										else: 
											m_all = factors(ns[0])

										while i_m < len(m_all):
											m = m_all[i_m]

											# initialize mut arrays for both 
											mutfound = np.zeros(n).reshape(1, n)

											mut1 = np.copy(mutfound)
											mut2 = np.copy(mutfound)

											# make up h_z values as well, pick just 1 value to start off 
											if h == 0.5:
												pop1_h = np.repeat(0.5, 1).reshape(1, 1)

											elif h == 8:
												pop1_h = np.random.uniform(low = 0.65, high = 1, size = 1).reshape(1, 1)

											elif h == 9:
												pop1_h = np.random.uniform(low = 0, high = 1, size = 1).reshape(1, 1)
											
											pop2_h = np.copy(pop1_h)

											# intitialize generation counter
											gen = 1
											WT_pheno1 = np.zeros(n).reshape(1, n)
											WT_pheno2 = np.copy(WT_pheno1)

											# evaluate how close the parents are to the optima
											# first, find the total_distance 
											origin = np.zeros(n)
											total_dist1 = np.linalg.norm(theta1 - origin)
											total_dist2 = np.linalg.norm(theta2 - origin)

											# find the remaning distance to reach the optima 
											remaining_dist1 = np.linalg.norm(WT_pheno1 - theta1)
											remaining_dist2 = np.linalg.norm(WT_pheno2 - theta2)

											# convert the remaining distance to percentage 
											percent1 = (remaining_dist1 / total_dist1) * 100
											percent2 = (remaining_dist2 / total_dist2) * 100

											gen = 1
											gen_max = 10000000
											# as long as one or both populations are adapting: 
											# ADAPTATION PHASE
											
											########################
											# ADAPTATION
											########################	

											while gen < gen_max:

												w1 = fitness(WT_pheno1, theta1, sigma_adapt)
												w2 = fitness(WT_pheno2, theta2, sigma_adapt)

												# generate a new mutations
												[hz1, newmut1] = mutate(alpha_adapt, n, h)
												[hz2, newmut2] = [hz1, newmut1]

												# calculate s
												s1 = compute_s(newmut1, WT_pheno1, theta1, sigma_adapt)
												s2 = compute_s(newmut2, WT_pheno2, theta2, sigma_adapt)

												# calculate h_fit
												hfit1 = compute_hfit(newmut1, WT_pheno1, theta1, sigma_adapt, s1, hz1)
												hfit2 = compute_hfit(newmut2, WT_pheno2, theta2, sigma_adapt, s2, hz2)

												# check if the muts fix
												POP1_did_it_fix = did_it_fix(s1, hfit1)
												POP2_did_it_fix = did_it_fix(s2, hfit2)
												
												# add them if they fixed
												[mut1, pop1_h] = evaluate_muts(POP1_did_it_fix, mut1, newmut1, pop1_h, hz1)
												[mut2, pop2_h] = evaluate_muts(POP2_did_it_fix, mut2, newmut2, pop2_h, hz2)

												# evaluate the new WT phenotype if the new mut fixes
												if POP1_did_it_fix == 1:
													WT_pheno1 = WT_pheno1 + newmut1
												if POP2_did_it_fix == 1:
													WT_pheno2 = WT_pheno2 + newmut2

												# find the remaning distance to reach the optima 
												remaining_dist1 = np.linalg.norm(WT_pheno1 - theta1)
												remaining_dist2 = np.linalg.norm(WT_pheno2 - theta2)

												# convert the remaining distance to percentage 
												percent1 = (remaining_dist1 / total_dist1) * 100
												percent2 = (remaining_dist2 / total_dist2) * 100

												gen += 1

												if percent1 <= circling_percentage and percent2 <= circling_percentage:
													print("Circling begins at gen", gen)
													break # end of while loop
											
											########################
											# CIRCLING
											########################									
											i = 0 
											gen_max = gen
											while i <= gen_max:

												w1 = fitness(WT_pheno1, theta1, sigma_adapt)
												w2 = fitness(WT_pheno2, theta2, sigma_adapt)

												# generate a new mutations
												[hz1, newmut1] = mutate(alpha_adapt, n, h)
												[hz2, newmut2] = [hz1, newmut1]

												# calculate s
												s1 = compute_s(newmut1, WT_pheno1, theta1, sigma_adapt)
												s2 = compute_s(newmut2, WT_pheno2, theta2, sigma_adapt)

												# calculate h_fit
												hfit1 = compute_hfit(newmut1, WT_pheno1, theta1, sigma_adapt, s1, hz1)
												hfit2 = compute_hfit(newmut2, WT_pheno2, theta2, sigma_adapt, s2, hz2)

												# check if the muts fix
												POP1_did_it_fix = did_it_fix(s1, hfit1)
												POP2_did_it_fix = did_it_fix(s2, hfit2)

												# add them if they fixed
												[mut1, pop1_h] = evaluate_muts(POP1_did_it_fix, mut1, newmut1, pop1_h, hz1)
												[mut2, pop2_h] = evaluate_muts(POP2_did_it_fix, mut2, newmut2, pop2_h, hz2)

												# evaluate the new WT phenotype if the new mut fixes
												if POP1_did_it_fix == 1:
													WT_pheno1 = WT_pheno1 + newmut1
												if POP2_did_it_fix == 1:
													WT_pheno2 = WT_pheno2 + newmut2

												gen += 1 # this is for the mutate function
												i += 1
											
											print("Circling ended at gen", i + gen_max - 1)
											final_gen = i + gen_max - 1
											print("FINAL_GEN", final_gen)

											#############################################
											# MODULARITY
											#############################################

											mut1 = add_modularity(mut1, modularity, m)
											mut2 = add_modularity(mut2, modularity, m)
											# print("m is", m)

											#############################################
											# CREATE THE ADAPTATION SUMMARY (of parents)
											#############################################
											
											one = np.ones(len(mut1)).reshape(len(mut1), 1)
											two = np.random.randint(low = 2, high = 3, size = len(mut2)).reshape(len(mut2), 1)
											pop1_or_pop2 = np.vstack((one, two))

											pop1_mut_n1 = mut1[:, 0].reshape(np.shape(mut1)[0], 1)
											pop2_mut_n1 = mut2[:, 0].reshape(np.shape(mut2)[0], 1)
											mut_n1_overall = np.round(np.vstack((pop1_mut_n1, pop2_mut_n1)), 3)
											
											pop1_mut_n2 = mut1[:, 1].reshape(np.shape(mut1)[0], 1)
											pop2_mut_n2 = mut2[:, 1].reshape(np.shape(mut2)[0], 1)
											mut_n2_overall = np.round(np.vstack((pop1_mut_n2, pop2_mut_n2)), 3)

											h_overall = np.vstack((pop1_h, pop2_h))

											# was modularity on or off
											modularity_table = np.repeat(modularity, np.shape(h_overall)[0]).reshape(np.shape(h_overall)[0], 1)

											# other stats like angle, alpha and dom 
											angle_table = np.repeat(round(angles[j]*180/math.pi,2), np.shape(h_overall)[0]).reshape(np.shape(h_overall)[0], 1)
											alpha_table = np.repeat(alpha_adapt, np.shape(h_overall)[0]).reshape(np.shape(h_overall)[0], 1)
											dom_table = np.repeat(h, np.shape(h_overall)[0]).reshape(np.shape(h_overall)[0], 1)

											#stack everyting in a table 
											adapt_summary_data = np.column_stack((angle_table, alpha_table, dom_table, pop1_or_pop2, mut_n1_overall, mut_n2_overall, h_overall, modularity_table))

											##################
											#HYBRIDIZATION / F1 - always choose parents from different populations 
											##################

											#loop over parallel evolution - on or off
											i_pevo = 0 

											if mut1[0,0] == 0 and mut1[0, 1] == 0:

												mut1 = np.delete(mut1, 0, axis = 0)
												mut2 = np.delete(mut2, 0, axis = 0)

												pop1_h = np.delete(pop1_h, 0, axis = 0)
												pop2_h = np.delete(pop2_h, 0, axis = 0)

												# remake the pop arrays here 
												pop1_chrom1 = np.ones(np.shape(mut1)[0])
												pop1_chrom2 = np.copy(pop1_chrom1)

												pop2_chrom1 = np.ones(np.shape(mut2)[0])
												pop2_chrom2 = np.copy(pop2_chrom1)


											while i_pevo < len(pevo):
												pevo_adapt = pevo[i_pevo]

												if pevo_adapt == 0: #pevo is off

													# hybrid genotype array is the same as hybrid_h
													hybrid_h = np.vstack((pop1_h, pop2_h))
													mut_hybrid = np.vstack((mut1, mut2))

													hybrid_all = hybrid_h

												else: #if pevo is on (there is parallel evolution), then do something different: 										

													#sum the mut matrices across columns 
													mut1sum = np.sum(mut1, axis = 1)
													mut2sum = np.sum(mut2, axis = 1)

													#find their intersection. this doesnt return indices, just returns the value(s) that is the same between the matrix and the intersection 
													intersect = np.intersect1d(mut1sum, mut2sum)
													
													#find indices in both mutsum arrays where mutsum has the same values as the intersect. columns that have those indices in the genotype matrices has the same mutations 
													pop1_pe_idx = np.flatnonzero(np.in1d(mut1sum, intersect))
													pop2_pe_idx = np.flatnonzero(np.in1d(mut2sum, intersect))

													################################
													# MUT ARRAYS
													################################

													mut1_pe = mut1[pop1_pe_idx]
													mut2_pe = mut2[pop2_pe_idx]

													mut1_del = np.delete(mut1, pop1_pe_idx, axis = 0)
													mut2_del = np.delete(mut2, pop2_pe_idx, axis = 0)
													
													mut_hybrid = np.vstack((mut1_pe, mut1_del, mut2_del))

													################################
													# H ARRAYS
													################################

													#find h values of parallel mutations
													pop1_h_pe = pop1_h[pop1_pe_idx]
													pop2_h_pe = pop2_h[pop2_pe_idx]

													#delete the h values of parallel muts, these are h values of non parallel loci 
													pop1_h_del = np.delete(pop1_h, pop1_pe_idx, axis = 0)
													pop2_h_del = np.delete(pop2_h, pop2_pe_idx, axis = 0)

													#make the hybrid h matrix. add the h values for 2nd pop 1st
													hybrid_h = np.vstack((pop1_h_pe, pop1_h_del, pop2_h_del))

													################################
													# hybrid genotypes
													################################

													parallel_loci = np.ones(np.shape(mut1_pe)[0])
													non_parallel_loci = np.append(pop1_h_del, pop2_h_del)

													hybrid_all = np.append(parallel_loci, non_parallel_loci)

												hybrid_pheno = np.dot(np.transpose(hybrid_all), mut_hybrid)

												# when h is variable, find the average h values 
												h_averaged = (np.mean(hybrid_h, axis = 0))

												# calculate the hybrid fitness in two optima and then assign them the higher of the two values 
												h_fit1 = np.array(fitness_F2(hybrid_pheno, theta1, sigma_adapt))
												h_fit2 = np.array(fitness_F2(hybrid_pheno, theta2, sigma_adapt))
												
												# make a random array as long as the h_fit1. We will replace the numbers here as we run the loop below. 
												h_fit = max(h_fit1, h_fit2)

												##################
												#PEVO METRIC  
												##################

												if pevo_adapt == 0:
													pevo_metric = 0
												else: 
													pevo_metric = 0.5 * (len(intersect) / int(np.shape(mut1)[0]) + len(intersect) / int(np.shape(mut2)[0])) * 100


												##################
												#TRAIT MISMATCH 
												##################
												mismatch = np.cross(WT_pheno2 - WT_pheno1,hybrid_pheno - WT_pheno1)/np.linalg.norm(WT_pheno2 - WT_pheno1)
												print("WT_pheno2 - WT_pheno1", WT_pheno2 - WT_pheno1)
												print("hybrid_pheno - WT_pheno1", hybrid_pheno - WT_pheno1)
												print("TOP", np.cross(WT_pheno2 - WT_pheno1,hybrid_pheno - WT_pheno1))

												##################
												#BACKCROSSES 
												##################
												
												# MAKE DIPLOID ONE VECTOR GENOTYPES 
												if pevo_adapt == 0: 
													[backcross_pheno1, backcross_fitness1] = make_backcrosses(mut1, mut2, None, None, None, pevo_adapt, "parent1", how_many_backcrosses, hybrid_h, mut_hybrid, theta1, theta2, sigma_adapt)
													[backcross_pheno2, backcross_fitness2] = make_backcrosses(mut2, mut1, None, None, None, pevo_adapt, "parent2", how_many_backcrosses, hybrid_h, mut_hybrid, theta1, theta2, sigma_adapt)

												else: 
													[backcross_pheno1, backcross_fitness1] = make_backcrosses(mut1, mut2, mut1_del, mut2_del, mut1_pe, pevo_adapt, "parent1", how_many_backcrosses, hybrid_h, mut_hybrid, theta1, theta2, sigma_adapt)
													[backcross_pheno2, backcross_fitness2] = make_backcrosses(mut2, mut1, mut1_del, mut2_del, mut1_pe, pevo_adapt, "parent2", how_many_backcrosses, hybrid_h, mut_hybrid, theta1, theta2, sigma_adapt)


												# find average backcross phenotype 
												backcross_pheno_mean1 = np.mean(np.array(backcross_pheno1), axis=0)
												backcross_pheno_mean2 = np.mean(np.array(backcross_pheno2), axis=0)

												# find average backcross fitness
												backcross_fit_mean1 = np.mean(np.array(backcross_fitness1), axis=0)
												backcross_fit_mean2 = np.mean(np.array(backcross_fitness2), axis=0)


												##################
												#HYBRIDIZATION / F2 - choose parents from the F1 
												##################

												if pevo == 1: 
													n_paralell = pop2_pe_idx.shape[0] # number of parallel loci  
													n_nonparallel = np.shape(mut1)[0] - n_paralell # number of non parallel loci 
												else: 
													n_paralell = 0 # number of parallel loci  
													n_nonparallel = np.shape(mut_hybrid)[0] - n_paralell # number of non parallel loci 

												# MAKE DIPLOID ONE VECTOR GENOTYPES 
												F2_list = list() # make an empty array to append later on 

												i = 0 
												while i < how_many_F2: 
													
													a = np.ones(n_paralell)
													b = np.random.choice([0, "h", 1], size = n_nonparallel, p = [0.25, 0.5, 0.25])
													F2 = np.append(a, b)
													F2_list.append(F2)
													
													i += 1
													
												# ADD ACTUAL H VALUES TO THE DIPLOID GENOTYPE	
												mut_hybrid_F2 = np.copy(mut_hybrid)
												hybrid_h = hybrid_h.flatten() # flatten into one dimension

												# add the actual h values into the F2 genotypes 
												for F2 in F2_list: 

													i = 0 
													while i < len(hybrid_h): 
														idx = np.where(F2 == "h")[0] # find the idx of "h"
														F2[idx] = hybrid_h[idx] # replace those with the real h values

														i += 1

												F2_list = [[float(x) for x  in F2] for F2 in F2_list]


												# CALCULATE PHENOTYPE
												pheno_list = []
												 
												for F2 in F2_list: 

													pheno = list(np.dot(F2, mut_hybrid))
													pheno_list.append(pheno)

												# find average F2 phenotype 
												F2_pheno_mean = np.mean(np.array(pheno_list), axis=0)


												# CALCULATE FITNESS
												fitness_list = []

												for pheno in pheno_list:

													w1 = np.array(fitness_F2(pheno, theta1, sigma_adapt))
													w2 = np.array(fitness_F2(pheno, theta2, sigma_adapt))

													# pick the higher one
													w = max(w1, w2)

													fitness_list.append(w)

												# find average F2 fitness
												F2_fit_mean = np.mean(np.array(fitness_list), axis=0)

												#######################
												# SAVE THE DATA 
												#######################

												# SAVE F2 AND BC VALUES FOR INDIVIDUALS, NOT MEANS
												individuals_table = np.hstack((
													np.repeat(rep+1, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(round(angles[j]*180/math.pi,2), how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(N_adapt, how_many_F2).reshape(how_many_F2, 1),
													np.repeat(sigma_adapt, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(u_adapt, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(alpha_adapt, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(opt_dist, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(h, how_many_F2).reshape(how_many_F2, 1),
													np.repeat(pevo_adapt, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(n, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(modularity, how_many_F2).reshape(how_many_F2, 1), 
													np.repeat(m, how_many_F2).reshape(how_many_F2, 1),
													np.repeat(h_fit, how_many_F2).reshape(how_many_F2, 1),
													np.array(fitness_list).reshape(how_many_F2, 1),
													np.array(backcross_fitness1).reshape(how_many_backcrosses, 1),
													np.array(backcross_fitness2).reshape(how_many_backcrosses, 1), 
													np.repeat(hybrid_pheno, how_many_F2).reshape(how_many_F2, 2),
													np.array(pheno_list).reshape(how_many_F2, 2),
													np.array(backcross_pheno1).reshape(how_many_backcrosses, 2),
													np.array(backcross_pheno2).reshape(how_many_backcrosses, 2)))										

												# SAVE FITNESS VALUES 
												fitness_table = np.hstack((rep+1, 
													round(angles[j]*180/math.pi,2),
													h_fit, 
													F2_fit_mean, 
													backcross_fit_mean1, 
													backcross_fit_mean2, 
													N_adapt, 
													sigma_adapt, 
													u_adapt, 
													alpha_adapt, 
													opt_dist, 
													h,
													np.mean(hybrid_h), # mean of all h values
													np.std(hybrid_h), # st dev of all h values
													pevo_adapt, 
													n, 
													modularity, 
													m)).reshape(1, 18)

												# SAVE PHENO VALUES 
												# compute mean and std of 
												param_list1 = [[rep + 1, round(angles[j]*180/math.pi,2)]]

												param_list2 = [[final_gen, N_adapt, 
														sigma_adapt, 
														u_adapt, 
														alpha_adapt, 
														opt_dist, 
														h, 
														pevo_adapt, 
														pevo_metric,
														n, 
														modularity, 
														m]]

												pheno_list = [WT_pheno1, 
													WT_pheno2,
													hybrid_pheno, 
													F2_pheno_mean, 
													backcross_pheno_mean1, 
													backcross_pheno_mean2,
													mismatch]
													

												pheno_list = [x.tolist() for x in pheno_list]
												
												# combine all lists 
												pheno_combined = param_list1 + pheno_list + param_list2

												import collections # to flatten the list 
												def flatten(x):
													if isinstance(x, collections.Iterable):
														return [a for i in x for a in flatten(i)]
													else:
														return [x]

												pheno_combined = flatten(pheno_combined)

												# here we calculate the number of cols in the phenotype array depending on the number of dimensions (n)
												ncol = len(pheno_list) * n + 13
												pheno_combined = np.asarray(pheno_combined).reshape(1, ncol)					
												
												# COMPUTE MEAN AND STD OF H VALUES
												hvals_combined = np.hstack((rep+1, 
													round(angles[j]*180/math.pi,2),
													final_gen, 
													N_adapt, 
													sigma_adapt, 
													u_adapt, 
													alpha_adapt, 
													opt_dist, 
													h,
													np.mean(pop1_h),
													np.std(pop1_h),
													np.mean(pop2_h), 
													np.std(pop2_h), 
													pevo_adapt, 
													pevo_metric,
													n, 
													modularity, 
													m)).reshape(1, 18)

												# save h mean and std values
												np.savetxt(fileHandle_A, hvals_combined, delimiter = ',')

												# save phenotype values
												np.savetxt(fileHandle_B, pheno_combined, delimiter = ',')

												# save adaptation summary
												adapt_summary = np.vstack((adapt_summary_data))
												np.savetxt(fileHandle_C, adapt_summary, fmt = '%.3f', delimiter = ',')

												# individual fitness and phenos of F2 and BC 
												np.savetxt(fileHandle_D, individuals_table, fmt = '%.3f', delimiter = ',')

												# save fitness values 
												np.savetxt(fileHandle_E, fitness_table, fmt = '%.3f', delimiter = ',')	

												# print an update
												if h == 9:
													print('N = %d, sigma = %.3f, u = %.3f, alpha = %.3f, opt_dist = %.2f, n=%d, angle=%r, rep=%d, h=%d, modularity=%0.1f, pevo=%0.1f, pevo_metric=%0.1f' %(N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n, round(angles[j]*180/math.pi,2), rep+1, h, modularity, pevo_adapt, pevo_metric)) 
												else:
													print('N = %d, sigma = %.3f, u = %.3f, alpha = %.3f, opt_dist = %.2f, n=%d, angle=%r, rep=%d, h=%s, modularity=%0.1f, pevo=%0.1f, pevo_metric=%0.1f' %(N_adapt, sigma_adapt, u_adapt, alpha_adapt, opt_dist, n, round(angles[j]*180/math.pi,2), rep+1, h, modularity, pevo_adapt, pevo_metric)) 

												#go to the next pevo value
												i_pevo += 1 

											# next factor for modularity 
											i_m += 1 

										# go to next rep
										rep += 1

									#go to next optimum
									j += 1

								#next dimension
								l += 1

							#next opt_dist
							i_opt_dist += 1 
							
							#next u_adapt 
							# i_u += 1 	
							
						#next alpha_adapt value
						i_alpha += 1

					# next modularity value
					i_modul += 1 

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

	# change to the final dir where we saved all the data frames  
	os.chdir(updated_path)

	# adaptation
	df_h = pd.read_csv("mutsummary.csv", header = None)
	df_h.columns = ['rep', 'angle', 'generation', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'paren1_mean_h', 'paren1_std_h', 'paren2_mean_h', 'paren2_std_h', 'pevo', 'pevo_metric', 'n', 'modularity', 'm'] 
	df_h.to_csv("mutsummary.csv")

	# adaptation
	df_adapt = pd.read_csv("adaptation.csv", header = None)
	df_adapt.columns = ['angle', 'alpha', 'dominance', 'pop1_or_pop2', 'mut_n1_overall', 'mut_n2_overall', 'h_overall', 'modularity'] 
	df_adapt.to_csv("adaptation.csv")

	if n == 2: 
		# do the renaming for phenotypes only if n = 2, for higher dims it gets complicated
		# renaming the phenodata is tricky since the number of cols will change depending on the number of dimensions 
		df_pheno = pd.read_csv("phenotypes.csv", header = None)
		df_pheno.columns = ['reps', 'angle', 'Parent1_1', 'parent1_2', 'parent2_1', 'parent2_2', 'F1_1', 'F1_2', 'F2_mean_1', 'F2_mean_2', 'backcross_pheno_mean1_1', 'backcross_pheno_mean1_2', 'backcross_pheno_mean2_1', 'backcross_pheno_mean2_2', 'mismatch', 'generation', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'pevo', 'pevo_metric', 'n', 'modularity', 'm'] # we wont change this first part
		df_pheno.to_csv("phenotypes.csv")

	elif n == 3: 
		# do the renaming for phenotypes only if n = 2, for higher dims it gets complicated
		# renaming the phenodata is tricky since the number of cols will change depending on the number of dimensions 
		df_pheno = pd.read_csv("phenotypes.csv", header = None)
		df_pheno.columns = ['reps', 'angle', 'Parent1_1', 'parent1_2', 'Parent1_3', 'parent2_1', 'parent2_2', 'parent2_3', 'F1_1', 'F1_2', 'F1_3',
		'F2_mean_1', 'F2_mean_2','F2_mean_3', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'pevo', 'pevo_metric', 'n', 'modularity', 'm'] # we wont change this first part
		df_pheno.to_csv("phenotypes.csv")

	elif n == 4: 
		df_pheno = pd.read_csv("phenotypes.csv", header = None)
		df_pheno.columns = ['reps', 'angle', 'Parent1_1', 'parent1_2', 'Parent1_3', 'parent1_4', 'parent2_1', 'parent2_2', 'parent2_3', 'parent2_4', 'F1_1', 'F1_2', 'F1_3', 'F1_4', 
		'F2_mean_1', 'F2_mean_2','F2_mean_3', 'F2_mean_4', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'pevo', 'pevo_metric', 'n', 'modularity', 'm'] # we wont change this first part
		df_pheno.to_csv("phenotypes.csv")

	elif n == 5: 
		df_pheno = pd.read_csv("phenotypes.csv", header = None)
		df_pheno.columns = ['reps', 'angle', 'Parent1_1', 'parent1_2', 'Parent1_3', 'parent1_4', 'parent1_5', 'parent2_1', 'parent2_2', 'parent2_3', 'parent2_4', 'parent2_4', 'F1_1', 'F1_2', 'F1_3', 'F1_4', 'F1_5', 
		'F2_mean_1', 'F2_mean_2','F2_mean_3', 'F2_mean_4', 'F2_mean_5', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'pevo', 'pevo_metric', 'n', 'modularity', 'm'] # we wont change this first part
		df_pheno.to_csv("phenotypes.csv")

	# fitMean
	df_fitMean = pd.read_csv("fitMean.csv", header = None)
	df_fitMean.columns = ['rep', 'angle', 'F1_fitness', 'F2_mean_fitness', 'backcross_fit_mean1', 'backcross_fit_mean1', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'mean_h', 'std_h', 'pevo', 'n', 'modularity', 'm']
	df_fitMean.to_csv("fitMean.csv")

	# F2 and BC
	df_fitMean = pd.read_csv("F2_BC.csv", header = None)
	df_fitMean.columns = ['rep', 'angle', 'N', 'sigma', 'u', 'alpha', 'opt_dist', 'h', 'pevo', 'n', 'modularity', 'm', 'F1_fitness', 'F2_fitness', 'BC1_fitness', 'BC2_fitness', 'F1_pheno1', 'F1_pheno2', 'F2_pheno1', 'F2_pheno2', 'BC_pheno1_1', 'BC_pheno1_2', 'BC_pheno2_1', 'BC_pheno2_2']
	df_fitMean.to_csv("F2_BC.csv")


	# add a little something at the end of the log file to show the run was successful
	f = open("success","w") # open a new text file

	f.close() # close the file, good practice! 

	print(updated_path)


######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    

start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
