# How dominant alleles influence progress toward speciation by natural selection
### A simulation study in python to understand speciation dynamics

Extrinsic post-zygotic isolating barriers result as a consequence of the alleles used for adaptation and occur when hybrids lack an appropriate ecological niche. Most models of this process consider only additive alleles, even though it is clear that alleles used for adaptation are not strictly additive. We used computer simulations based on Fisher's geometric model of adaptation to investigate how the availability of dominant alleles for adaptation in parental populations affects the phenotype and fitness of their hybrids. We find that dominant alleles typically reduce hybrid fitness under parallel selection because hybrids tend to have transgressive phenotypes that are maladaptive in the common parental niche. Under divergent selection, dominance typically increases hybrid fitness because chance events cause unequal dominance among populations, causing hybrids to resemble one parent more than the other. Although in organisms with many traits dominance reduces hybrid fitness. Our results indicate that dominant alleles facilitate progress toward speciation via parallel natural selection, and might either increase or decrease the efficacy of divergent selection for speciation.

Simulations from this study are now published in the journal Evolution: https://academic.oup.com/evolut/article/76/12/2846/6975846

---------------------------------------------------------------------

The `maintext_community.py` script runs the simulations and relies on helper functions defined in `utilities.py`. Upon execution, it generates several output files, each of which is described below. All output files include informative column headers, detailing both the simulated parameters and the resulting metrics.

1. `phenotypes.csv`: Contains the phenotypes of the parents and the both classes of hybrids.
2. `adaptation.csv`: This file keeps track of the mutations that arise during adaptation. Based on it, it is possible to track each mutation down to the population it arose and its effect size. Each row is for a single mutation that arose and fixed in the population. Lost mutations aren't tracked.
3. `mutsummary.csv`: This file gives more information about the mutation events such as the mean dominance coefficients and standard deviations of mutations in the parental populations. 
4. `fitMean.csv`: This file contains information about the fitness of the parents as well as both classes of hybrids. Note that we consider the mean fitness for the `BC` and `F2` populations. 
5. `F2_BC.csv`: This file has information about the fitness and phenotypes of individuals in the `BC` and `F2` populations. Note that depending on how many individuals you choose to include, this file may get quite large.
