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
