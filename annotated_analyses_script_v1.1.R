#
#  v.1 		November 7, 2014
#  v.1.1 	February 8, 2014
#
#  This file contains all annotated R scripts to repeat our analyses. It demonstrates the 
#  implementation of our analyses, and contains comments on results.  However, running the
#  analyses verbatim from this script (in serie) takes very long.  Therefore, we executed 
#  the analyses that use the function corHMM::corDISC() on the OSCAR computer cluster at 
#  the Centre for Computation and Visualization at Brown University. To facilitate 
#  following all our steps, we provide 11 checkpoints: saves of intermediate results at the 
#  end of each section that can be loaded to start the next section of the analyses.
# 
#
#  All analyses implemented and annotated by Jurriaan M. de Vos.  Feel free to contact 
#  with questions or remarks at jurriaan_devos@brown.edu
#
#
#  The script is divided into 4 parts, corresponding with four objectives:
#   - objective 1.  Ascertain whether Zanne et al. analyzed vessel cross-sectional area or diameter
#   - objective 2.  Replicate the main results published in the main text of the paper
#   - objective 3.  Investigate the effects of climate grooming on pathway results
#   - objective 4.  Explore behavior of the "% pathway" approach using simulations
#  Each part is subdivided in numbered Tasks.

#############################################################################################

#  get ready
	rm(list=ls())
	
	# set work directory
	basedir <- "/foo/bar/github/"  ## change this to where you copied the files from github
	# as an example, I copied them here:
	basedir <- "~/Desktop/pers_work/project_zanne/github/"
	setwd(basedir)
	
	# load functions
	library(corHMM)
	library(stringr) # just for function str_replace_all()
	library(phytools) # just for simulations 
	source('R_source_custom_functions.R')
	# this source file contains the following functions:
	# pickSimul() is a function to simply compare two corDISC result objects based on $AIC 
	# 	element to print which one is better fitting
	# getPathways() is a function that takes as input a 4x4 transition rate matrix and 
	# 	calculates from it the % contribution of the different pathways from (0,0) to (1,1),
	# 	as clarified to us upon contact with Jeremy Beaulieu and Brian O'Meara around 
	# 	May 19, 2014
	# getPathways_print() is a function that takes as input a corDISC result object and 
	# 	uses it to print the % pathway and other results (such as trait lability and 
	# 	persistence times)
	# countFirstTransitions() is a function to count from a simmap trees the number of 
	# 	times the state (1,1) (i.e. final state) was reached from states (0,0), (0,1), and
	# 	(1,0). It only counts the first transitions (looking from the root forward in 
	# 	time), not the total number of transitions to (1,1).
	# fileListToResultTable() is a function to take a list of files that contain corDISC 
	# 	results, load them one by one, and summarize them in a table which is written to a 
	# 	file
	# resultTableToSelectedModels() is a function that reads a table created by function 
	# 	fileListToResultTable() and returns similar table file but with only best models 
	# 	for each analysis type. 



##########################################################################################################
### objective 1.  Ascertain whether Zanne et al. analyzed vessel cross-sectional area rather than diameter
##########################################################################################################

	#########################
	#  task 1.1  read and format data and tree objects for analysis
	#########################

# read data
	clim   <- read.csv("./dryad/MinimumFreezingExposure.csv", stringsAsFactors=FALSE)
	vessel <- read.table("./dryad/GlobalVessel_A_data.txt", header=TRUE)
	tree   <- read.tree("./dryad/Vascular_Plants_rooted.dated.tre")
	woody  <- read.csv("./dryad/GlobalWoodinessDatabase.csv")

# adjust and cross-check data
	
	# format climate data
	clim$species <- str_replace_all(clim$species, " ", "_")
	clim$Freeze.tmin.lo[clim$Freeze.tmin.lo=="FreezingExposed"] <- 1
	clim$Freeze.tmin.lo[clim$Freeze.tmin.lo=="FreezingUnexposed"] <- 0
	clim$climate <- as.numeric(clim$Freeze.tmin.lo)
	clim <- clim[clim$species %in% tree$tip.label, ]
	clim <- data.frame(species=clim$species, climate=clim$climate, row.names=clim$species)

	# throw out non-woody species from vessel data and format
	woody$species <- str_replace_all(woody$gs, " ", "_")
	nonwoodyspecies <- woody$species[woody$woodiness == "H"]
	nonwoodyvesselspecies <- vessel$Binomial[vessel$Binomial %in% nonwoodyspecies]
	vessel <- vessel[!vessel$Binomial %in% nonwoodyvesselspecies, ]
	vessel <- vessel[vessel$Binomial %in% tree$tip.label, ]
	vessel <- data.frame(species=vessel$Binomial, A=vessel$A, D=2*sqrt(vessel$A/pi), row.names=vessel$Binomial)
		# Note, vessel$A is the surface area data; vessel$D is these values calculated into diameters

	# drop tips from tree for which there is no data
	tipstodrop <- tree$tip.label[!(tree$tip.label %in% vessel$species  &  tree$tip.label %in% clim$species)]
	tree <- drop.tip(tree, tip=tipstodrop)


# make input data objects for corDISC analysis

	data_D <- data.frame(species=tree$tip.label, climate=NA, vessel_D=NA)
	for(spec in data_D$species) {  # ok looping is silly but it works
		data_D$climate[data_D$species == spec] <- clim$climate[clim$species == spec]
		data_D$vessel_D[data_D$species == spec] <- vessel$D[vessel$species == spec]
	}
	small <- data_D$vessel_D < 0.044
	data_D$vessel_D[small] <- 1
	data_D$vessel_D[!small] <- 0
	
	data_A <- data.frame(species=tree$tip.label, climate=NA, vessel_A=NA)
	for(spec in data_A$species) { 
		data_A$climate[data_A$species == spec] <- clim$climate[clim$species == spec]
		data_A$vessel_A[data_A$species == spec] <- vessel$A[vessel$species == spec]
	}
	small <- data_A$vessel_A < 0.044
	data_A$vessel_A[small] <- 1
	data_A$vessel_A[!small] <- 0

	# make an index matrix to allow simultaneous transitions	
	rate_mat <- matrix(data=c(NA, 1, 2, 9,  3, NA, 10, 4,  5, 11, NA, 6,  12, 7, 8, NA), byrow=FALSE, nrow=4)
	rownames(rate_mat) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
	colnames(rate_mat) <- rownames(rate_mat)

	# ----------------------
	# --  Checkpoint 01 
	# ----------------------
	save(file="./intermediate_results/01_ObjectsToAnalyzeVesselAvsD.Rdata", list=c("tree", "data_D", "data_A", "rate_mat"))


	#########################
	#  task 1.2  execute analyses
	#########################

# run corDISC analyses on A vs D data

	rm(list=ls())
	load(file="./intermediate_results/01_ObjectsToAnalyzeVesselAvsD.Rdata")
	corDISCresult_null_A <- corDISC(phy=tree, data=data_A, rate.mat=NULL, model="ARD")
	corDISCresult_null_D <- corDISC(phy=tree, data=data_D, rate.mat=NULL, model="ARD")
	corDISCresult_ratemat_A <- corDISC(phy=tree, data=data_A, rate.mat=rate_mat, model="ARD")
	corDISCresult_ratemat_A <- corDISC(phy=tree, data=data_D, rate.mat=rate_mat, model="ARD")
	# ----------------------
	# --  Checkpoint 02 
	# ----------------------
	save(file="./intermediate_results/01_ObjectsToAnalyzeVesselAvsD.Rdata", list=c("tree", "data_D", "data_A", "rate_mat"))
	save.image(file="./intermediate_results/02_results_image_AvsD.Rdata")


	#########################
	#  task 1.3  summarize results
	#########################


# analyze results I: determine whether simultaneous transitions are supported based on AIC
	rm(list=ls())
	load(file="./intermediate_results/02_results_image_AvsD.Rdata")
	source('R_source_custom_functions.R')

	pickSimul(nos=corDISCresult_null_A, s=corDISCresult_ratemat_A)
	# gives:
	# Simultaneous changes not supported; use model null

	pickSimul(nos=corDISCresult_null_D, s=corDISCresult_ratemat_D)
	# gives:
	# Simultaneous changes not supported; use model null


# analyze results II: get pathways, state lability, and persistence times

	getPathways_print(corDISCresult_null_A)
	 #  returns:
	 #
	 #  Summary of corDISC result:
	 #  --------------------------
	 #  percent trait-first  : 88.76
	 #  percent simultaneous : 0
	 #  percent climate-first: 11.24
	 #  
	 #  climate- vs trait-lability: 1
	 #  
	 #  persistence time clim.ances & trait.ances: 10.011
	 #  persistence time clim.ances & trait.deriv: 41.215
	 #  persistence time clim.deriv & trait.ances: 0.005
	 #  persistence time clim.deriv & trait.deriv: 34.92
	 #  
	 #  number taxa with clim.ances & trait.ances: 17
	 #  number taxa with clim.ances & trait.deriv: 450
	 #  number taxa with clim.deriv & trait.ances: 0
	 #  number taxa with clim.deriv & trait.deriv: 369
	 #  
	 #  fraction of taxa with freezing-adapted trait state: 0.98

	getPathways_print(corDISCresult_null_D)
	 #  returns:
	 #  Summary of corDISC result:
	 #  --------------------------
	 #  percent trait-first  : 55.1
	 #  percent simultaneous : 0
	 #  percent climate-first: 44.9
	 #  
	 #  climate- vs trait-lability: 5.248
	 #  
	 #  persistence time clim.ances & trait.ances: 52.075
	 #  persistence time clim.ances & trait.deriv: 5.244
	 #  persistence time clim.deriv & trait.ances: 21.378
	 #  persistence time clim.deriv & trait.deriv: 11.812
	 #  
	 #  number taxa with clim.ances & trait.ances: 393
	 #  number taxa with clim.ances & trait.deriv: 74
	 #  number taxa with clim.deriv & trait.ances: 169
	 #  number taxa with clim.deriv & trait.deriv: 200
	 #  
	 #  fraction of taxa with freezing-adapted trait state: 0.328

	### Conclusion: the analysis using cross sectional area (A) gives results very 
	#   close to the numbers Zanne et al.'s [we: (88.76, 0, 11.24) vs. Zanne et al.: (82.7, 0, 17.3)]
	#   The numbers are slightly different either way, the cause for which is unknown to us.
	#   We do note that the number of taxa are also slightly different (836 in our analyses based on 
	#   Zanne et al.'s data files vs. 860 reported in their paper)


##########################################################################
### objective 2.  Replicate the main results published in the main text of the paper
##########################################################################
	rm(list=ls())
	basedir <- "/foo/bar/github/"  ## change this to where you copied the files from github
	setwd(basedir)


	#########################
	#  task 2.1  read and format data and tree objects
	#########################

# read data
	clim   <- read.csv("./dryad/MinimumFreezingExposure.csv", stringsAsFactors=FALSE)
	tree   <- read.tree("./dryad/Vascular_Plants_rooted.dated.tre")
	woody  <- read.csv("./dryad/GlobalWoodinessDatabase.csv")
	phenol <- read.csv("./dryad/GlobalLeafPhenologyDatabase.csv")

# adjust and cross-check data
	# format species in Genus_species format in all datasets
	phenol$Spec <- str_replace_all(phenol$Binomial, " ", "_")
	phenol$Spec <- factor(phenol$Spec)
	woody$Spec <- str_replace_all(woody$gs, " ", "_")
	woody$Spec <- factor(woody$Spec)

	# format climate data
	clim$Spec <- str_replace_all(clim$species, " ", "_")
	clim$Freeze.tmin.lo[clim$Freeze.tmin.lo=="FreezingExposed"] <- 1
	clim$Freeze.tmin.lo[clim$Freeze.tmin.lo=="FreezingUnexposed"] <- 0
	clim$climate <- as.numeric(clim$Freeze.tmin.lo)
	clim <- clim[clim$Spec %in% tree$tip.label, ]
	clim <- data.frame(Spec=clim$Spec, climate=clim$climate, row.names=clim$Spec)

	# woodiness data
	woody <- woody[!woody$woodiness=="variable", ]  # throw out variable species
	woody <- woody[woody$Spec %in% tree$tip.label,] # throw out species not in the tree
	woody$dat <- 0
	woody$dat[woody$woodiness == "H"] <- 1  # score 1 for harbaceous taxa

	# throw out species not in the tree nor belonging in the data; add scoring
	# leaf phenology	
	phenol <- phenol[phenol$Phenology!="D_EV",]  # throw out variable species (p. 11 supplementary material).
	phenol <- phenol[phenol$Spec %in% tree$tip.label,] # throw out species not in the tree
	phenol <- phenol[! phenol$Spec %in% woody$Spec[woody$woodiness == "H"], ] # throw out non-woody speces
											# keep species of which it is not known that they are herbaceous
	phenol$dat <- 0							# score traits: 1 if deciduous
	phenol$dat[phenol$Phenology == "D"] <- 1

	# tree

	# extract angiosperms and subtrees (rosids, asterids, magnolids, monocots, "the rest", from the full tree)
	node <- which(tree$node.label == "Angiospermae")+length(tree$tip.label)
	tree <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Superrosidae")+length(tree$tip.label)
	rosids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Superasteridae")+length(tree$tip.label)
	asterids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Magnoliidae")+length(tree$tip.label)
	magnolids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Monocotyledoneae")+length(tree$tip.label)
	monocots <- extract.clade(tree, node=node)
	tipsToDrop <- c(rosids$tip.label, asterids$tip.label, magnolids$tip.label, monocots$tip.label)
	resttree <- drop.tip(tree, tipsToDrop)

	# make an index matrix to allow simultaneous transitions	
	rate_mat <- matrix(data=c(NA, 1, 2, 9,  3, NA, 10, 4,  5, 11, NA, 6,  12, 7, 8, NA), byrow=FALSE, nrow=4)
	rownames(rate_mat) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
	colnames(rate_mat) <- rownames(rate_mat)

	# make data objects

	# phenol
		includedSpecies <- as.character(phenol$Spec[ (phenol$Spec %in% tree$tip.label) & (phenol$Spec %in% clim$Spec) ])
		includedSpecies <- includedSpecies[ !includedSpecies %in% woody$Spec[woody$dat==1] ]  # making sure to not include herbaceous species
		phenol_data <- data.frame(
			spec=includedSpecies, 
			climate=clim[ clim$Spec %in% includedSpecies, "climate"], 
			phenol=phenol[ phenol$Spec %in% includedSpecies, "dat"],
			row.names=includedSpecies
		)
	# woody
		includedSpecies <- as.character(woody$Spec[ (woody$Spec %in% tree$tip.label) & (woody$Spec %in% clim$Spec) ])
		woody_data <- data.frame(
			spec=includedSpecies, 
			climate=clim[ clim$Spec %in% includedSpecies, "climate"], 
			woody=woody[ woody$Spec %in% includedSpecies, "dat"],
			row.names=includedSpecies
		)

	# ----------------------
	# --  Checkpoint 03 
	# ----------------------

	save(list=c("tree", "asterids", "rosids", "magnolids", "monocots", "resttree", "woody_data", "phenol_data", "rate_mat"), file="./intermediate_results/03_ObjectsToAnalyzePhenolWoody.Rdata")

	#########################
	#  task 2.2  run the analyses
	#########################

	#  run analyses  (we ran these on a computer cluster)
		rm(list=ls())
		load(file="./intermediate_results/03_ObjectsToAnalyzePhenolWoody.Rdata")

		# Woodiness analyses analyses
		phy <- rosids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_rosids_woody.Rdata")
		phy <- asterids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_asterids_woody.Rdata")
		phy <- monocots
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_monocots_woody.Rdata")
		phy <- magnolids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_magnolids_woody.Rdata")
		phy <- resttree
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_resttree_woody.Rdata")

		phy <- rosids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_rosids_woody.Rdata")
		phy <- asterids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_asterids_woody.Rdata")
		phy <- monocots
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_monocots_woody.Rdata")
		phy <- magnolids
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_magnolids_woody.Rdata")
		phy <- resttree
		data <- woody_data[ woody_data[,1] %in% phy$tip.label, ]
		phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% woody_data[,1] ] )	
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_resttree_woody.Rdata")

		# Leaf phenology analyses
		phy <- tree
		data <- phenol_data[ phenol_data[,1] %in% phy$tip.label, ]
		result <- corDISC(phy=phy, data=data, rate.mat=NULL, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_magnolids_woody.Rdata")
		result <- corDISC(phy=phy, data=data, rate.mat=rate_mat, model="ARD")
		save(result,  file="./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_resttree_woody.Rdata")

	# ----------------------
	# --  Checkpoint 04 
	# ----------------------
	# all individual analyses were saved to file
	

	#########################
	#  task 2.3  summarize results for leaf phenology
	#########################


	load("./intermediate_results/04_FilesResultsPhenolWoody/corDISC_NULL_repeat_phenol.Rdata")
	repeat_NULL_phenol <- result
	load("./intermediate_results/04_FilesResultsPhenolWoody/corDISC_rate_mat_repeat_phenol.Rdata")
	repeat_ratemat_phenol <- result

	source('R_source_custom_functions.R')

	pickSimul(nos=repeat_NULL_phenol, s=repeat_ratemat_phenol)
	# gives:
	# Simultaneous changes  supported; use model ratemat

	getPathways_print(repeat_ratemat_phenol)
	 # Summary of corDISC result:
	 # --------------------------
	 # percent trait-first  : 39.28
	 # percent simultaneous : 17.64
	 # percent climate-first: 43.07
	 # 
	 # climate- vs trait-lability: 0.689
	 # 
	 # persistence time clim.ances & trait.ances: 49.84
	 # persistence time clim.ances & trait.deriv: 11.71
	 # persistence time clim.deriv & trait.ances: 17.752
	 # persistence time clim.deriv & trait.deriv: 44.065
	 # 
	 # number taxa with clim.ances & trait.ances: 1524
	 # number taxa with clim.ances & trait.deriv: 174
	 # number taxa with clim.deriv & trait.ances: 424
	 # number taxa with clim.deriv & trait.deriv: 471
	 # 
	 # fraction of taxa with freezing-adapted trait state: 0.249
	 # 
	 # 
	### in the paper it was: 
	 # 	Trait first 36.7%
	 #  Simultaneous 14.6%
	 #  Climate first 48.7%
	# so our results are very similar (but, again, not identical)

	#########################
	#  task 2.4.  summarize results for woodiness
	#########################

	# get the results as objects
	fileList <- list.files("./intermediate_results/04_FilesResultsPhenolWoody/")
	fileList <- fileList[grepl("_woody", fileList)]
	for (name in fileList) {
		load(paste("./intermediate_results/04_FilesResultsPhenolWoody/", name, sep=""))
		assign(strsplit(name, split=".Rdata")[[1]][1], result)
	}

	# model selection: are simultaneous changes supported?
	pickSimul(corDISC_NULL_repeat_asterids_woody, corDISC_rate_mat_repeat_asterids_woody)
		# Simultaneous changes  supported; use model ratemat
	pickSimul(corDISC_NULL_repeat_rosids_woody, corDISC_rate_mat_repeat_rosids_woody)
		# Simultaneous changes  supported; use model ratemat
	pickSimul(corDISC_NULL_repeat_magnolids_woody, corDISC_rate_mat_repeat_magnolids_woody)
		# Simultaneous changes not supported; use model null
	pickSimul(corDISC_NULL_repeat_monocots_woody, corDISC_rate_mat_repeat_monocots_woody)
		# Simultaneous changes not supported; use model null
	pickSimul(corDISC_NULL_repeat_resttree_woody, corDISC_rate_mat_repeat_resttree_woody)
		# Simultaneous changes not supported; use model null

	# the best 5-partition model would have 12+12+8+8+8=48 parameters.  The paper states a 40 parameter model was choosen.
	# Perhaps model selection was performed jointly, not allowing for simultaneous changes in only part of the tree.
	# Let's use the joint-likelihood expression and compare calculated AIC scores
	lik_NULL <- sum(corDISC_NULL_repeat_asterids_woody$loglik, corDISC_NULL_repeat_rosids_woody$loglik, corDISC_NULL_repeat_magnolids_woody$loglik, corDISC_NULL_repeat_monocots_woody$loglik, corDISC_NULL_repeat_resttree_woody$loglik)  
	lik_ratemat <- sum(corDISC_rate_mat_repeat_asterids_woody$loglik, corDISC_rate_mat_repeat_rosids_woody$loglik,  corDISC_rate_mat_repeat_magnolids_woody$loglik,  corDISC_rate_mat_repeat_monocots_woody$loglik, corDISC_rate_mat_repeat_resttree_woody$loglik)
	AIC_NULL    <- 2*5*8 - 2*lik_NULL
	AIC_ratemat <- 2*5*12 - 2*lik_ratemat
	dAIC <- abs(AIC_NULL-AIC_ratemat)
	if(AIC_NULL < AIC_ratemat) {cat("best is NULL; dAIC =", dAIC)} else {cat("best is ratematall; dAIC =", dAIC)}
	# gives: best is NULL; dAIC = 5.654045
	# Thus, this way there is indeed support preference for no-simultaneous-transitions
	# (note that the optimal model changes when we adjust climate scoring below).

	# to get to the "the relative likelihood of the different pathways out of the woody and freezing-unexposed state and into the herbaceous and freezing-exposed state", there are two possible meanings of "The weighted average (by clade diversity)": 
	# -1- calculating the weighted averages of the transition rates and get the pathways from them
	# -2- calculating 5 pathways, and averaging those directly. 
	# We try both approaches to see which one yields results as stated in the paper.
	# Note that the latter approach weights the FIRST transitions into the freezing by the EXTANT number of species. Hypothetically, if a group diversifies rapidly after transitioning for the first time into the freezing zone, that transition would be given more importance than a transition in a group that also colonized the freezing zone but did not radiate rapidly upon arrival.  Approach -1- seems more in line with the spirit of the paper, but the results in the paper are replicated using approach -2- (see below), so we assume that that is how they were calculated.


	# -1- : calculate pathways via calculating weighted average of transition rates: 
	# This leads to results that are somewhat different from Zanne et al.
		rateproduct <- 
			corDISC_NULL_repeat_asterids_woody$solution * dim(corDISC_NULL_repeat_asterids_woody$data)[1] +
			corDISC_NULL_repeat_rosids_woody$solution * dim(corDISC_NULL_repeat_rosids_woody$data)[1] +
			corDISC_NULL_repeat_magnolids_woody$solution * dim(corDISC_NULL_repeat_magnolids_woody$data)[1] +
			corDISC_NULL_repeat_monocots_woody$solution * dim(corDISC_NULL_repeat_monocots_woody$data)[1] 
			#		+	corDISC_NULL_repeat_resttree_woody$solution * dim(corDISC_NULL_repeat_resttree_woody$data)[1]
		tips <-  
			dim(corDISC_NULL_repeat_asterids_woody$data)[1] +
			dim(corDISC_NULL_repeat_rosids_woody$data)[1] +
			dim(corDISC_NULL_repeat_magnolids_woody$data)[1] +
			dim(corDISC_NULL_repeat_monocots_woody$data)[1] 
			# +
			# dim(corDISC_NULL_repeat_resttree_woody$data)[1]
		pars <- rateproduct/tips

		# get pathways
		pathways <- getPathways(pars)
		traitfirst <- pathways[1]
		simult <- pathways[2]
		climfirst <- pathways[3]
		
		# lability ratio:  
		climratesum  <- sum(pars[1,3], pars[2,4], pars[1,4], pars[3,1], pars[4,2], pars[4,1], na.rm=TRUE)
		traitratesum <- sum(pars[1,2], pars[3,4], pars[1,4], pars[2,1], pars[4,3], pars[4,1], na.rm=TRUE)
		clim_vs_trait_labil <- climratesum / traitratesum
	
		# persistence times
		persis00  <- 1/rowSums(pars, na.rm=TRUE)[1]
		persis01  <- 1/rowSums(pars, na.rm=TRUE)[2]
		persis10  <- 1/rowSums(pars, na.rm=TRUE)[3]
		persis11  <- 1/rowSums(pars, na.rm=TRUE)[4]
	
		phrase <- paste(
			"Summary of corDISC result:",
			"\n--------------------------",
			"\npercent trait-first  : ", round(traitfirst, 5),
			"\npercent simultaneous : ", round(simult, 2),
			"\npercent climate-first: ", round(climfirst, 2),
			"\n",
			"\nclimate- vs trait-lability: ", round(clim_vs_trait_labil, 3),
			"\n",
			"\npersistence time clim.ances & trait.ances: ", round(persis00, 3), 
			"\npersistence time clim.ances & trait.deriv: ", round(persis01, 3),
			"\npersistence time clim.deriv & trait.ances: ", round(persis10, 3), 
			"\npersistence time clim.deriv & trait.deriv: ", round(persis11, 3),
			"\n",
			sep=""
		)
		cat(phrase)
	# gives:
	#  Summary of corDISC result when going via weighted average :
	#  --------------------------
	#  percent trait-first  : 66.33
	#  percent simultaneous : 0
	#  percent climate-first: 33.67
	#  
	#  climate- vs trait-lability: 6.501
	#  
	#  persistence time clim.ances & trait.ances: 75.413
	#  persistence time clim.ances & trait.deriv: 4.225
	#  persistence time clim.deriv & trait.ances: 24.246
	#  persistence time clim.deriv & trait.deriv: 14.235

	# These results that are quite a bit different from Zanne et al., who reported 58% and 42% from trait-first and climate-first, respectively.


	# -2- : calculate pathways via directly averaging the results for the individual trees: 
	# This leads to results that are very much in line with the ones reported by Zanne et al.
		pathways_all <- NULL
		specs <- NULL
	for(result in c("corDISC_NULL_repeat_asterids_woody", "corDISC_NULL_repeat_rosids_woody", "corDISC_NULL_repeat_magnolids_woody", "corDISC_NULL_repeat_monocots_woody", "corDISC_NULL_repeat_resttree_woody")) 
	{
		pars <- get(result)$solution
		pathways <- getPathways(pars)
		pathways_all <- rbind(pathways_all, pathways)
		specs <- c(specs, dim(get(result)$data)[1])
	}
		colnames(pathways_all) <- c("trait_first", "simultaneous", "climate_first")
		rownames(pathways_all) <- c("asterids", "rosids", "magnolids", "monocots", "resttree")
		pathways_all
			#             trait_first simultaneous climate_first
			#  asterids  8.228868e+01            0      17.71132
			#  rosids    7.497492e+01            0      25.02508
			#  magnolids 0.000000e+00            0     100.00000
			#  monocots  0.000000e+00            0     100.00000
			#  resttree  1.122390e-10            0       0.00000
			## note the "resttree" result does not sum to 100% because 
			##  00->01 is 10^18 times smaller than 01->00 and 
			##  and while 00->10 occurs, 10->11 = 0). Thus, species
			##  hardly ever reach 11 from 00...  Thus, the resttree must have been excluded:
		pathways_all <- pathways_all[1:4,]
		specs <- specs[1:4]
		traitfirst <- sum(pathways_all[,1]*(specs/sum(specs)))
		simult <- sum(pathways_all[,2]*(specs/sum(specs)))
		climfirst <- sum(pathways_all[,3]*(specs/sum(specs)))
		traitfirst		# 56.4113 %
		simult			#  0   %
		climfirst		# 43.5887 %
		#  this is much closer to 58%, 0%, and 42%, respectively, in the paper, which we thus assume is the approach followed.



##########################################################################
### objective 3.  Investigate the effects of climate grooming on pathway results
##########################################################################
	rm(list=ls())
	basedir <- "/foo/bar/github/"  ## change this to where you copied the files from github
	setwd(basedir)

	#########################
	#  task 3.1  read and format data and tree objects
	#########################
	
# read data

	clim <- read.table("./handling_climate_data/climate_data_grooming_levels_JdV.txt", header=TRUE)  
	# this is the climate data we obtained, across different grooming levels and percentile thresholds

	vessel <- read.table("./dryad/GlobalVessel_A_data.txt", header=TRUE)
	woody  <- read.csv("./dryad/GlobalWoodinessDatabase.csv")
	phenol <- read.csv("./dryad/GlobalLeafPhenologyDatabase.csv")
	tree   <- read.tree("./dryad/Vascular_Plants_rooted.dated.tre")
	# these are the datasets from Dryad.

# adjust and cross-check data
	# format species in Genus_species format in all datasets
	phenol$Spec <- str_replace_all(phenol$Binomial, " ", "_")
	phenol$Spec <- factor(phenol$Spec)
	woody$Spec <- str_replace_all(woody$gs, " ", "_")
	woody$Spec <- factor(woody$Spec)
	vessel$Spec <- factor(vessel$Binomial)
	clim$Spec <- clim$Spec


	# throw out species not in the tree nor belonging in the data; add scoring

	# leaf phenology	
	phenol <- phenol[phenol$Phenology!="D_EV",]  # throw out variable species (p. 11 supplementary material).
	phenol <- phenol[phenol$Spec %in% tree$tip.label,] # throw out species not in the tree
	phenol <- phenol[! phenol$Spec %in% woody$Spec[woody$woodiness == "H"], ] # throw out non-woody speces
											# keep species of which it is not known that they are herbaceous
	phenol$dat <- 0							# score traits: 1 if deciduous
	phenol$dat[phenol$Phenology == "D"] <- 1


	# vessel data
	vessel <- vessel[! vessel$Spec %in% woody$Spec[woody$woodiness == "H"], ] # throw out non-woody speces
	vessel <- vessel[vessel$Binomial %in% tree$tip.label, ] # throw out species not in the tree
	vessel <- data.frame(Spec=vessel$Binomial, D=2*sqrt(vessel$A/pi), row.names=vessel$Binomial)
	vessel$dat <- 0
	vessel$dat[vessel$D < 0.044] <- 1  # score 1 for small vessel using 0.044 mm threshold for diameter
		# Note,  vessel$D is surface area data calculated into diameters


	# woodiness data
	woody <- woody[!woody$woodiness=="variable", ]  # throw out variable species
	woody <- woody[woody$Spec %in% tree$tip.label,] # throw out species not in the tree
	woody$dat <- 0
	woody$dat[woody$woodiness == "H"] <- 1  # score 1 for harbaceous taxa

	# throwout climate records of species that do not matter
	clim <- clim[(clim$Spec %in% tree$tip.label &
		(clim$Spec %in% vessel$Spec |
		 clim$Spec %in% woody$Spec |
		 clim$Spec %in% phenol$Spec)), ]  # keep species that are in the tree and in either of the three data sets
	# summarize the percentage of species with few unique localities
	(sum(clim[,"Records_groomed"]<10) / dim(clim)[1])*100
	# gives 12.6182 % of the species that are part of an analysis have 4 or less unique localities

	# give the list of species that change scoring at percentile threshold 2.5% due to grooming:
	# species scored as tropical by Zanne et al. and temperate by us:
	as.character(clim$Spec[clim$Tmin_all_perc0025 > 0 & clim$Tmin_groomed_perc0025 < 0])
	# species scored as temperate by Zanne et al. and tropical by us:
	as.character(clim$Spec[clim$Tmin_all_perc0025 < 0 & clim$Tmin_groomed_perc0025 > 0])
	
	
	# extract angiosperms from tree
	node <- which(tree$node.label == "Angiospermae")+length(tree$tip.label)
	tree <- extract.clade(tree, node=node)
	# drop tips from tree for which there is no data
	tipstodrop <- tree$tip.label[ !(tree$tip.label %in% clim$Spec &
		(tree$tip.label %in% vessel$Spec |
		 tree$tip.label %in% woody$Spec |
		 tree$tip.label %in% phenol$Spec))]
	tree <- drop.tip(tree, tip=tipstodrop)
		
	# summarize for how many species the threshold percentile matters
	overlap <- clim[clim$Spec %in% tree$tip.label, ]
	total <- dim(overlap)[1]  # number of species in tree and climate data: 13536   12882
	doubtful <- sum(overlap$Tmin_all_percmin <= 0   &  overlap$Tmin_all_percmax > 0)  # 6452
	round((doubtful/total)*100,2) 
	#  gives 50.09
	# thus, 50% of the species in the tree for which data is available occur both in freezing and non-freezing conditions

	# transform climate data: score 0 for when Tmin ≤ 0; score 1 when Tmin >0
	TminColumns <- colnames(clim)[grepl("Tmin", colnames(clim))]  # for the columns that are climate data
	for(i in TminColumns){
		notFreezing <- clim[, i] > 0
		freezing <- clim[, i] <= 0
		clim[notFreezing, i] <- 0		
		clim[freezing, i] <- 1
	}

	#########################
	#  task 3.2  create the objects as input for corDISC analyses
	#########################
	# make 126 input data objects, by combining 42 climate data sets with the woody, phenol and vessel data. 
	#	42 climate data sets:  7 climate thresholds (at minimum Tmin, 2.5 percentile (as Zanne et al.), 10 percentile, 50 percentile (median), 90 percentile, 97.5 percentile, maximum), for two grooming levels (all records; all records after excluding duplicates, imprecise and doubtful coordinates), for three minimal number of records (1, 3, 10).

# make conduit size ("vessel") objects
	for (percentile in c("percmin", "perc0025", "perc0100", "percmed", "perc0900", "perc0975", "percmax")) {
		for (nmin in c(1, 3, 10)) {
			for(groomingLevel in c("all", "groomed")) {
				# to assist in indexing:
				includedSpecies <- as.character(vessel$Spec[ (vessel$Spec %in% tree$tip.label) & (vessel$Spec %in% clim$Spec[clim$Records_all >= nmin ]) ])
				includedSpecies <- includedSpecies[ !includedSpecies %in% woody$Spec[woody$dat==1] ]  # making sure to not include herbaceous species
				climateDataName <- paste("Tmin", groomingLevel, percentile, sep="_")
		
				assign(paste("data", groomingLevel, nmin, percentile, "vessel", sep="_"), 
					data.frame(
						spec=includedSpecies, 
						percentile=clim[ clim$Spec %in% includedSpecies, climateDataName], 
						vessel=vessel[ vessel$Spec %in% includedSpecies, "dat"],
						row.names=includedSpecies
					)
				)
			}
		}
	}
# make leaf phenology ("phenol") objects
	for (percentile in c("percmin", "perc0025", "perc0100", "percmed", "perc0900", "perc0975", "percmax")) {
		for (nmin in c(1, 3, 10)) {
			for(groomingLevel in c("all", "groomed")) {
				# to assist in indexing:
				includedSpecies <- as.character(phenol$Spec[ (phenol$Spec %in% tree$tip.label) & (phenol$Spec %in% clim$Spec[clim$Records_all >= nmin ]) ])
				includedSpecies <- includedSpecies[ !includedSpecies %in% woody$Spec[woody$dat==1] ]  # making sure to not include herbaceous species
				climateDataName <- paste("Tmin", groomingLevel, percentile, sep="_")
		
				assign(paste("data", groomingLevel, nmin, percentile, "phenol", sep="_"), 
					data.frame(
						spec=includedSpecies, 
						percentile=clim[ clim$Spec %in% includedSpecies, climateDataName], 
						phenol=phenol[ phenol$Spec %in% includedSpecies, "dat"],
						row.names=includedSpecies
					)
				)
			}
		}
	}

# make leaf woodiness ("woody") objects
	for (percentile in c("percmin", "perc0025", "perc0100", "percmed", "perc0900", "perc0975", "percmax")) {
		for (nmin in c(1, 3, 10)) {
			for(groomingLevel in c("all", "groomed")) {
				# to assist in indexing:
				includedSpecies <- as.character(woody$Spec[ (woody$Spec %in% tree$tip.label) & (woody$Spec %in% clim$Spec[clim$Records_all >= nmin ]) ])
				climateDataName <- paste("Tmin", groomingLevel, percentile, sep="_")
		
				assign(paste("data", groomingLevel, nmin, percentile, "woody", sep="_"), 
					data.frame(
						spec=includedSpecies, 
						percentile=clim[ clim$Spec %in% includedSpecies, climateDataName], 
						woody=woody[ woody$Spec %in% includedSpecies, "dat"],
						row.names=includedSpecies
					)
				)
			}
		}
	}

# make subtrees on which the the woodiness data can be analyzed

	node <- which(tree$node.label=="Superrosidae")+length(tree$tip.label)
	rosids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Superasteridae")+length(tree$tip.label)
	asterids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Magnoliidae")+length(tree$tip.label)
	magnolids <- extract.clade(tree, node=node)
	node <- which(tree$node.label=="Monocotyledoneae")+length(tree$tip.label)
	monocots <- extract.clade(tree, node=node)
	tipsToDrop <- c(rosids$tip.label, asterids$tip.label, magnolids$tip.label, monocots$tip.label)
	resttree <- drop.tip(tree, tipsToDrop)

# make the index matrix to analyze simultaneous rates

	rate_mat <- matrix(data=c(NA, 1, 2, 9,  3, NA, 10, 4,  5, 11, NA, 6,  12, 7, 8, NA), byrow=FALSE, nrow=4)
	rownames(rate_mat) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
	colnames(rate_mat) <- rownames(rate_mat)

#  save objects needed for analysis: the trees and the data and rate_mat

	# ----------------------
	# --  Checkpoint 05 / 06 
	# ----------------------
	
	save(file="./intermediate_results/05_ObjectsDataClimateScoringSensitivity.Rdata", list = ls()[grepl("data_", ls())])
	save(file="./intermediate_results/06_ObjectsTreerunnerClimateScoringSensitivity.Rdata", list=c("tree", "rosids", "asterids", "magnolids", "monocots", "resttree", "rate_mat"))


	#########################
	#  task 3.3  run analyses
	#########################
	rm(list=ls())
	load("./intermediate_results/05_ObjectsDataClimateScoringSensitivity.Rdata")
	inputlist <- ls()
	load("./intermediate_results/06_ObjectsTreerunnerClimateScoringSensitivity.Rdata")

	##  each data object created above is analyzed using corHMM::corDISC()
	##  the result of each analysis is stored in its own result file
	##  running the analyses this way takes long; you may want to split up the loops 
	##	and run them in parallel on a cluster
	##  
	##  the set of analyess using the index matrix rate_mat specifying that simultaneous changes are allowed
	for(i in 1:length(inputlist)) {
		# get data objects
		data <- get(inputlist[i])
		# select rate mat
		rate.mat <- rate_mat

		# vessel and phenology data are analyzed across the whole tree
		if(grepl("_vessel", inputlist[i]) | grepl("_phenol", inputlist[i])) { 
			# match tree to data
			phy <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data[,1]])
			# run analysis
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD")
			# store result object with meaningful file name
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_tree_", inputlist[i], ".Rdata", sep=""))
		}	
		
		# woody data are analyzed across each of the 5 subtrees
		if(grepl("_woody", inputlist[i])) { 

			phy <- drop.tip(asterids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_asterids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(rosids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_rosids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(monocots, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_monocots_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(magnolids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_magnolids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(resttree, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_ratematall_resttree_", inputlist[i], ".Rdata", sep="")) #save analysis
		}
	}


	##  the set of analyess using no index matrix so simultaneous changes are not allowed
	for(i in 1:length(inputlist)) {
		# get data objects
		data <- get(inputlist[i])
		# select rate mat
		rate.mat <- NULL

		# vessel and phenology data are analyzed across the whole tree
		if(grepl("_vessel", inputlist[i]) | grepl("_phenol", inputlist[i])) { 
			# match tree to data
			phy <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data[,1]])
			# run analysis
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD")
			# store result object with meaningful file name
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_tree_", inputlist[i], ".Rdata", sep=""))
		}	
		
		# woody data are analyzed across each of the 5 subtrees
		if(grepl("_woody", inputlist[i])) { 

			phy <- drop.tip(asterids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_asterids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(rosids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_rosids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(monocots, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_monocots_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(magnolids, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_magnolids_", inputlist[i], ".Rdata", sep="")) #save analysis

			phy <- drop.tip(resttree, tree$tip.label[!tree$tip.label %in% data[,1]])  # match tree to data
			data <- data[data[,1]%in%phy$tip.label,]  # match data to tree
			result <- corDISC(phy=phy, data=data, rate.mat=rate.mat, model="ARD") # run analysis
			save(result, file=paste("./intermediate_results/07_FilesResultsClimateScoringSensitivity/
corHMM_result_NULL_resttree_", inputlist[i], ".Rdata", sep="")) #save analysis
		}
	}
	# ----------------------
	# --  Checkpoint 07 
	# ----------------------
	# all individual analyses were saved to file


	#########################
	#  task 3.4  compile results in a table
	#########################

	rm(list=ls())
	source('R_source_custom_functions.R')

	##  calculate each individual result object into a line in a table containing all relevant aspects of the results
	fileList <- list.files("./intermediate_results/07_FilesResultsClimateScoringSensitivity", full.names=TRUE)
	fileListToResultTable(fileList, outputfile="./intermediate_results/08_outputTable.txt")

	# for each analysis in the table, grab the relevant lines and do model selection and averaging
	res <- read.table("./intermediate_results/08_outputTable.txt", header=TRUE, stringsAsFactors=FALSE)

	resultTableToSelectedModels(res, subset=TRUE, trait="vessel", outputfile="results_ClimateScoringSensitivity.txt")
	resultTableToSelectedModels(res, subset=TRUE, trait="phenol", outputfile="results_ClimateScoringSensitivity.txt", append=TRUE)
	resultTableToSelectedModels(res, subset=TRUE, trait="woody", treetype="compound", testtype="direct", outputfile="results_ClimateScoringSensitivity.txt", append=TRUE)

	# ----------------------
	# --  Checkpoint 08
	# ----------------------
	# all results were saved to a table
	# that table is called "results_ClimateScoringSensitivity.txt"

	#########################
	#  task 3.5.  summarize results into a figure 1A for the main text
	#########################
	rm(list=ls())
	source('R_source_custom_functions.R')
	res <- read.table("results_ClimateScoringSensitivity.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

	pdf("./figures/Figure_ClimateScoringSensitivity.pdf", width=8.9/2.54, height=8.9/2.54, fonts="sans")

##  plot settings		
	simultcol <- "#555555"
	polcol <- "white"
	bgcol <- "grey85"
	simulcex <- 0.5
	mtextline <- -0.85
	phenol_col <- "green"
	vessel_col <- "royalblue1"
	vessel_col <- "blue"
	woody_col <- "red"
	cexdata <- 0.75

## set-up plotting area
	par(mar=c(1.5,1.5,0.5, 0.5))
	plot(NULL, bty="o",
		xlim=c(0,100), ylim=c(0,100), 
		type="n", xaxt="n", yaxt="n")
	axis(side=1, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=-3, tck=-0.02)
	axis(side=2, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=2, tck=-0.02)
	## add axes text 
		# in the paper, it is stated "contribution of different evolutionary solutions"
	mtext(side=1, line=0.5, "Inferred % climate-first pathway", cex=0.6)
	mtext(side=2, line=0.8, "Inferred % trait-first pathway", cex=0.6)

## draw polygons for areas in which one trait wins
	#background
	xx <- c(  0, 0, 100)
	yy <- c(0, 100, 0)
	polygon(x=xx, y=yy, col=bgcol, border= bgcol, lwd=2)

	xx <- c(  0, 44, 29.333, 0)
	yy <- c(100, 56, 41.333, 56)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_trait <- max(xx)-(max(xx)-min(xx))/2
	ytext_trait <- max(yy)-(max(yy)-min(yy))/2
	
	xx <- c(56, 100, 56, 41.333)
	yy <- c(44,   0,  0, 29.333)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_clim <- max(xx)-(max(xx)-min(xx))/2
	ytext_clim <- max(yy)-(max(yy)-min(yy))/2

	xx <- c(0, 44, 29.333,  0)
	yy <- c(0,  0, 29.333, 44)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_simul <- max(xx)-(max(xx)-min(xx))/2
	ytext_simul <- max(yy)-(max(yy)-min(yy))/2

	# add text
	textcex <- 0.66
	text(x=3, y=66, "Trait-first", cex= textcex, pos=4)
	text(x=51, y=18, "Climate-first", cex= textcex, pos=4)
	text(x=3, y=18, "Simultaneous", cex= textcex, pos=4)
	text(x=18, y=37, "Gray zone", cex= textcex, pos=4)


# select data to plot
	res <- res[res$groom=="all" & res$minclimdat==1, ]

# plot data
	# format data to transform into colour
	res$perc <- factor(res$perc, levels=c("percmin", "perc0025", "perc0100", "percmed", "perc0900", "perc0975", "percmax"))

	maxalph <- 255
	pch <- 21

	phenolcol <- vector()
	woodycol <- vector()
	vesselcol <- vector()

	# greens plus making all colors:
	data <- res[res$trait=="phenol", ]
		for(n in 1:dim(data)[1]) {
		i <- 255
		j <- 255 - ((maxalph / length(levels(res$perc)))*n)
		phenolcol[n] <- rgb(red=0, green=i, blue=0, alpha=j, maxColorValue=255)
		woodycol[n] <- rgb(red=i, green=0, blue=0, alpha=j, maxColorValue=255)
		vesselcol[n] <- rgb(red=0, green=0, blue=i, alpha=j, maxColorValue=255)
	}

	nudge <- 0.5  # to "jitter" points fully overlaying in the 100% trait-first category

	which <- which(data$traitfirst > 99.9)
	data$traitfirst[which] <- data$traitfirst[which] + nudge
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col="white", bg="white", cex= cexdata))
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=phenol_col, bg=phenolcol, cex= cexdata))

	# reds:
	data <- res[res$trait=="woody", ]
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col="white", bg="white", cex= cexdata))
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=woody_col, bg=woodycol, cex= cexdata))

	# blues:
	data <- res[res$trait=="vessel", ]
	which <- which(data$traitfirst > 99.9 & data$perc=="perc0900")
	data$climfirst[which] <- data$climfirst[which] - nudge
	which <- which(data$traitfirst > 99.9 & data$perc=="perc0975")
	data$climfirst[which] <- data$climfirst[which] + nudge
	which <- which(data$traitfirst > 96 & data$perc=="percmax")
	data$climfirst[which] <- data$climfirst[which] + nudge*1.5

	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col="white", bg="white", cex= cexdata))
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=vessel_col, bg=vesselcol, cex= cexdata))

## add zanne result (z) and corrigendum (z*)
	cexdot <- 1.2

	col <- phenol_col
	points(x=48.7, y=36.7, cex=cexdot, pch="+", col=col)

	col <- woody_col
	points(x=42.0, y=58.0, cex=cexdot, pch="+", col=col)

	col <- vessel_col
	points(x=17.3, y=82.7, cex=cexdot, pch="+", col=col)
	text(x=44.9, y=55.1, cex=cexdot, label="+", col=col)	
	text(x=44.9, y=55.1, cex=cexdot, label="  *", col=col)	

# add legend
	base <- 45
	textcex <- 0.6
	text(x=base, y=100, "+", cex=0.666, pos=4, font=1)
	text(x=base+3.9, y=100, ": Point estimates in Zanne et al.", cex=textcex, pos=4)
	text(x=base, y=95, "+*", cex=0.666, pos=4, font=1)
	text(x=base+3.9, y=95, ": Point estimate in corrigendum", cex=textcex, pos=4)

	x <- c(base+seq(4,18,14/6))
	y <- rep(90, length(x))
	points(x, y, pch=pch, cex=0.75, col=phenol_col, bg=phenolcol)
	text(x=base+16, y=mean(y), ": Deciduous leaves", cex=textcex, pos=4)	
	y <- y-5
	points(x, y, pch=pch, cex=0.75, col=vessel_col, bg=vesselcol)
	text(x=base+16, y=mean(y), ": Narrow conduits", cex=textcex, pos=4)	
	y <- y-5
	points(x, y, pch=pch, cex=0.75, col=woody_col, bg=woodycol)
	text(x=base+16, y=mean(y), ": Herbaceous growth form", cex=textcex, pos=4)	
	y <- y-3	
	text(x=x, y=y, labels=c("min", "2.5", "10", "50", "90", "97.5", "max"), cex=0.4, srt=270, adj=0)
	text(x=base+16, y=mean(y-2), ": Percentile threshold", cex=textcex, pos=4)	

# Add "A"
	mtext("A", side=1, line=0, adj=-2/30)

dev.off()
	

	#########################
	#  task 3.6.  summarize all results into 6 figures [part of the suplemental]
	#########################
	# In each of six panels, a figure is plotted similar to Figure 1A (i.e., with 8 
	# percentile thresholds to score "freezing-exposed" for each of the data types).  
	# Each is based on the results of analyses that used different procedures to cleanse 
	# the GBIF data. Either all records were included, or GBIF records were groomed 
	# (removing duplicate records, those with low geographic precision, with locality 
	# identical to a herbarium or centroid of a major political region - see the annotated 
	# python script for details); and either no minimum sample size was used, or species 
	# needed at least 3 or 10 unique localities to be included.  
	
	# It is not straight forward to decide what the best cleansing method and percentile 
	# thresholds are, but it is evident from this figure that these decisions have a 
	# tremendous effect on the results.

	rm(list=ls())
	source('R_source_custom_functions.R')
	res <- read.table("results_ClimateScoringSensitivity.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

	pdf("./figures/Figure_ClimateScoringSensitivity_all.pdf", width=(8.9/2.54)*2, height=(8.9/2.54)*3, fonts="sans")

## set-up plotting area
	par(mar=c(1.5,1.5,0.5, 0.5), mfrow=c(3,2), cex=1)

	for(panel in 1:6) {
		
##  plot settings		
	simultcol <- "#555555"
	polcol <- "white"
	bgcol <- "grey85"
	simulcex <- 0.5
	mtextline <- -0.85
	phenol_col <- "green"
	vessel_col <- "royalblue1"
	vessel_col <- "blue"
	woody_col <- "red"
	cexdata <- 0.75

	plot(NULL, bty="o",
		xlim=c(0,100), ylim=c(0,100), 
		type="n", xaxt="n", yaxt="n")
	axis(side=1, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=-3, tck=-0.02)
	axis(side=2, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=2, tck=-0.02)
	## add axes text 
		# in the paper, it is stated "contribution of different evolutionary solutions"
	mtext(side=1, line=0.5, "Inferred % climate-first pathway", cex=0.6)
	mtext(side=2, line=0.8, "Inferred % trait-first pathway", cex=0.6)

## draw polygons for areas in which one trait wins
	#background
	xx <- c(  0, 0, 100)
	yy <- c(0, 100, 0)
	polygon(x=xx, y=yy, col=bgcol, border= bgcol, lwd=2)

	xx <- c(  0, 44, 29.333, 0)
	yy <- c(100, 56, 41.333, 56)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_trait <- max(xx)-(max(xx)-min(xx))/2
	ytext_trait <- max(yy)-(max(yy)-min(yy))/2
	
	xx <- c(56, 100, 56, 41.333)
	yy <- c(44,   0,  0, 29.333)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_clim <- max(xx)-(max(xx)-min(xx))/2
	ytext_clim <- max(yy)-(max(yy)-min(yy))/2

	xx <- c(0, 44, 29.333,  0)
	yy <- c(0,  0, 29.333, 44)
	polygon(x=xx, y=yy, col= polcol, border=NA)
	xtext_simul <- max(xx)-(max(xx)-min(xx))/2
	ytext_simul <- max(yy)-(max(yy)-min(yy))/2

	# add text
	textcex <- 0.66
	text(x=3, y=66, "Trait-first", cex= textcex, pos=4)
	text(x=51, y=18, "Climate-first", cex= textcex, pos=4)
	text(x=3, y=18, "Simultaneous", cex= textcex, pos=4)
	text(x=18, y=37, "Gray zone", cex= textcex, pos=4)


# select data to plot

	if(panel == 1) {res2 <- res[res$groom=="all" & res$minclimdat==1, ]}
	if(panel == 2) {res2 <- res[res$groom=="groomed" & res$minclimdat==1, ]}
	if(panel == 3) {res2 <- res[res$groom=="all" & res$minclimdat==3, ]}
	if(panel == 4) {res2 <- res[res$groom=="groomed" & res$minclimdat==3, ]}
	if(panel == 5) {res2 <- res[res$groom=="all" & res$minclimdat==10, ]}
	if(panel == 6) {res2 <- res[res$groom=="groomed" & res$minclimdat==10, ]}
	

## add zanne result (z) and corrigendum (z*)
	cexdot <- 1.2

	col <- phenol_col
	points(x=48.7, y=36.7, cex=cexdot, pch="+", col=col)

	col <- woody_col
	points(x=42.0, y=58.0, cex=cexdot, pch="+", col=col)

	col <- vessel_col
	points(x=17.3, y=82.7, cex=cexdot, pch="+", col=col)
	text(x=44.9, y=55.1, cex=cexdot, label="+", col=col)	
	text(x=44.9, y=55.1, cex=cexdot, label="  *", col=col)	


# plot data
	# format data to transform into colour
	res2$perc <- factor(res2$perc, levels=c("percmin", "perc0025", "perc0100", "percmed", "perc0900", "perc0975", "percmax"))
	pch <- c(15, 17, 19, 0, 2, 1, 4)

	cexdata <- cexdata-0.1
	
	# greens:
	data <- res2[res2$trait=="phenol", ]
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=phenol_col, cex=cexdata))

	# reds:
	data <- res2[res2$trait=="woody", ]
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=woody_col, cex= cexdata))

	# blues:
	data <- res2[res2$trait=="vessel", ]
	with(data, points(traitfirst[order(perc)] ~ climfirst[order(perc)], pch=pch, col=vessel_col, cex= cexdata))


	# add legend
	base <- 45
	textcex <- 0.6
	ystart=90
	if(panel == 1) {
		text(x=base, y= ystart, "+", cex=0.666, pos=4, font=1)
		text(x=base+3.9, y= ystart, ": Point estimates in Zanne et al.", cex=textcex, pos=4)
		text(x=base, y= ystart-5, "+*", cex=0.666, pos=4, font=1)
		text(x=base+3.9, y= ystart-5, ": Point estimate in corrigendum", cex=textcex, pos=4)

		x <- c(base+seq(4,18,14/6))
		y <- rep(ystart-10, length(x))
		points(x, y, pch=pch, cex=0.75, col=phenol_col)
		text(x=base+16, y=mean(y), ": Deciduous leaves", cex=textcex, pos=4)	
		y <- y-5
		points(x, y, pch=pch, cex=0.75, col=vessel_col)
		text(x=base+16, y=mean(y), ": Narrow conduits", cex=textcex, pos=4)	
		y <- y-5
		points(x, y, pch=pch, cex=0.75, col=woody_col)
		text(x=base+16, y=mean(y), ": Herbaceous growth form", cex=textcex, pos=4)	
		y <- y-3	
		text(x=x, y=y, labels=c("min", "2.5", "10", "50", "90", "97.5", "max"), cex=0.4, srt=270, adj=0)
		text(x=base+16, y=mean(y-2), ": Percentile threshold", cex=textcex, pos=4)	
	}
	if(panel %in% c(1, 3, 5)) {
		text(x=base, y=ystart+10, "Raw GBiF records", cex=textcex, pos=4, font=2)
	}
	if(panel %in% c(2, 4, 6)) {
		text(x=base, y=ystart+10, "Cleansed GBiF records", cex=textcex, pos=4, font=2)
	}
	
	if(panel %in% c(1, 2)) {
		text(x=base, y=ystart+5, "Minimum 1 record per species", cex=textcex, pos=4, font=2)
	}
	if(panel %in% c(3, 4)) {
		text(x=base, y=95, "Minimum 3 records per species", cex=textcex, pos=4, font=2)
	}
	if(panel %in% c(5, 6)) {
		text(x=base, y=95, "Minimum 10 records per species", cex=textcex, pos=4, font=2)
	}

}
dev.off()
	






##########################################################################
### objective 4.  Explore behavior of the "% pathway" approach using simulations
##########################################################################

	rm(list=ls())
	basedir <- "/foo/bar/github/"  ## change this to where you copied the files from github
	setwd(basedir)
	source("R_source_custom_functions.R")

	#########################
	#  task 4.1  generate simmap simulations
	#########################

	# objective of the simulations is to address
	#  -  if there is a strong preference for a particular pathway in the underlying model (e.g. trait-first), is it correctly inferred?
	#  -  if pathways are equally likely, are pathways correctly inferred to be equally likely?
	#
	# For interpretation, we guide us by the statement by Zanne et al.: "[leaf phenology] was also far more likely to evolve as a response to a change in environment rather than arising before encountering freezing (that is, climate occupancy evolved first)", which was made in reference to a 12% higher preference of the climate-first pathway relative to the trait-first pathway (36.7% trait-first vs. 48.7% climate-first).

	# as tree to simulate along, we take the tree used for the vessel analyses.
	clim   <- read.csv("./dryad/MinimumFreezingExposure.csv", stringsAsFactors=FALSE)
	clim$species <- str_replace_all(clim$species, " ", "_")
	vessel <- read.table("./dryad/GlobalVessel_A_data.txt", header=TRUE)
	tree   <- read.tree("./dryad/Vascular_Plants_rooted.dated.tre")
	tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% vessel$Binomial & tree$tip.label %in% clim$species)])
	# tree has 836 species

	## we simulate without simultaneous transitions. (When all 12 transitions are allowed and all rates are equal, the expectation is that simulataneous transitions are more common because non-simultaneous changes require two steps).

	# We simululate under two scenarios:
	#  -1- all rates equal
	#  -2- 3x preference for trait-first trajectory
		
	# SCENARIO 1:  q matrix with all rates equal, no simultaneous transitions
	q <- matrix(nrow=4, byrow=TRUE, data=c(
			0   , 0.01, 0.01, 0   , 
			0.01, 0   , 0   , 0.01,
			0.01, 0   , 0   , 0.01, 
			0   , 0.01, 0.01, 0   ))
	diag(q) <- -rowSums(q)
	# these trees contain no (stochasitcially) prefered trajectory. rate 0.01 was chosen for simulation because such low-ish rates yield more "first" transitions, because under very high rates, the transitions will already occur on the branches coming from the root, hence, there will only be 2 first-transitions, and we don't give the method a fair chance.
	#  we only keep trees in which the actual number of first transitions to state (1,1) is within 2.5% of what we expect 
	#    (namely 50% trait-first under scenario 1, and 75% trait-first under scenario 2)

	i <- 0
	n <- 0
	sims <- list()
	# generate trees
	while(n<200){
		i <- i+1
		cat(paste("\n", i, "...  "))
		x <- sim.history(tree, q, anc="1", nsim=1)
		paths <- countFirstTransitions(x)	
		paths <- (paths[1:3]/paths[4]) * 100
		if (abs(paths[1] - paths[3]) < 5) {
			cat(paste("we have a winner in round", i))
			n <- n+1
			sims[[n]] <- x
			cat(paste("!!!; sims total is", n))
		} else {cat(paste("difference", abs(paths[1] - paths[3]), "%")) }
	}
	sims_equal <- sims

	# SCENARIO 2:  q matrix with 75% preferred trait-first pathway (3x more likely), no simultaneous transitions
	q <- matrix(nrow=4, byrow=TRUE, data=c(
			#00	  #01   #10   #11				# read these as from column-state to row-state
			0    , 1/300, 0.01  , 0   , # 00 
			0.01 , 0    , 0    , 0.01,	# 01
			0.01 , 0    , 0    , 0.01, 	# 10
			0    , 0.01 , 1/300, 0   ))	# 11 
	diag(q) <- -rowSums(q)
	i <- 0
	n <- 0
	sims <- list()
	# generate trees
	while(n<200){
		i <- i+1
		cat(paste("\n", i, "...  "))
		x <- sim.history(tree, q, anc="1", nsim=1)
		paths <- countFirstTransitions(x)	
		paths <- (paths[1:3]/paths[4]) * 100
		if (paths[1] < 77.5 & paths[1] > 72.5) {
			cat(paste("we have a winner in round", i))
			n <- n+1
			sims[[n]] <- x
			cat(paste("!!!; sims total is", n))
		} else {cat(paste("difference", abs(paths[1] - paths[3]), "%")) }
	}
	sims_traitfirst <- sims

	# ----------------------
	# --  Checkpoint 09
	# ----------------------
	save(file="./intermediate_results/09_simulatedTrees.Rdata", list=c("sims_traitfirst", "sims_equal"))  


	#########################
	#  task 4.2  run corDISC analyses
	#########################
	rm(list=ls())
	load("./intermediate_results/09_simulatedTrees.Rdata")
	# run analyses  (this is slow, may want to roll out loops and run on a cluster)
	for (i in 1:length(sims_equal)) {
			# format data
			x <- sims_traitfirst[[i]]
	  		simdat <- data.frame(spec=x$tip.label, clim=NA, trait=NA)
	  		simdat[,"clim"][x$states%in%c(1,2)] <- "0"
	  		simdat[,"clim"][x$states%in%c(3,4)] <- "1"
	  		simdat[,"trait"][x$states%in%c(1,3)] <- "0"
	  		simdat[,"trait"][x$states%in%c(2,4)] <- "1"
	  		# run
	  		result <- corDISC(x, simdat, rate.mat=NULL, model="ARD")
	  		# save
	  		save(result, file=paste("./intermediate_results/10_FilesResultsSimulations/simResults_equal_", formatC(i, width=3, format="d", flag="0"), ".Rdata", sep=""))
	}

	for (i in 1:length(sims_traitfirst)) {
			# format data
			x <- sims_traitfirst[[i]]
	  		simdat <- data.frame(spec=x$tip.label, clim=NA, trait=NA)
	  		simdat[,"clim"][x$states%in%c(1,2)] <- "0"
	  		simdat[,"clim"][x$states%in%c(3,4)] <- "1"
	  		simdat[,"trait"][x$states%in%c(1,3)] <- "0"
	  		simdat[,"trait"][x$states%in%c(2,4)] <- "1"
	  		# run
	  		result <- corDISC(x, simdat, rate.mat=NULL, model="ARD")
	  		# save
	  		save(result, file=paste("./intermediate_results/10_FilesResultsSimulations/simResults_traitfirst_", formatC(i, width=3, format="d", flag="0"), ".Rdata", sep=""))
	}
	# ----------------------
	# --  Checkpoint 10
	# ----------------------
	# all individual results saved to file

		
	###################################################
	#  task 4.3  process analyses and store in objects pathways_inf
	###################################################

	rm(list=ls())

	load("./intermediate_results/09_simulatedTrees.Rdata")
	source("R_source_custom_functions.R")

	# get results from equal-rates trajectory (no simultaneous changes)
	# first inferred from inferred rate matrix, then from true counts of first transitions in the simmap tree
	fileList <- list.files("./intermediate_results/10_FilesResultsSimulations", full.names=TRUE)
	fileList <- fileList[grepl("simResults_equal", fileList)]
	pathways_inf <- matrix(data=NA, nrow=0, ncol=3)
	colnames(pathways_inf) <- c("traitfirst", "simultaneous", "climatefirst")
	for (i in 1:length(fileList)) {
		# load result
		load(fileList[i])

		# get pathway
		pars <- result$solution
		pathways <- getPathways(pars)

		# store pathway
		pathways_inf <- rbind(pathways_inf, pathways)
	}

	pathways_true <- matrix(data=NA, nrow=0, ncol=3)
	colnames(pathways_true) <- c("traitfirst", "simultaneous", "climatefirst")
	# count truenumber of first transitions into % pathways
	for(i in 1:length(sims_equal)) {
		cat("\nworking on tree", i)
		x <- sims_equal[[i]]

		result <- countFirstTransitions(x)
		pathways <- (result[1:3]/result[4]) * 100  # [1:3] are the number of times a pathway is taken;  [4] is total number of first transitions
		pathways_true <- rbind(pathways_true, pathways)
	}


	# get results from biased trajectory (no simultaneous changes)
	# first inferred from inferred rate matrix, then from true counts of first transitions in the simmap tree
par(mfrow=c(1,2))
	fileList <- list.files("./intermediate_results/10_FilesResultsSimulations", full.names=TRUE)
	fileList <- fileList[grepl("simResults_traitfirst", fileList)]
	pathways_biased_inf <- matrix(data=NA, nrow=0, ncol=3)
	colnames(pathways_biased_inf) <- c("traitfirst", "simultaneous", "climatefirst")
	for (i in 1:length(fileList)) {
		# load result
		load(fileList[i])

		# get pathway
		pars <- result$solution
		pathways <- getPathways(pars)

		# store pathway
		pathways_biased_inf <- rbind(pathways_biased_inf, pathways)
	}

	pathways_biased_true <- matrix(data=NA, nrow=0, ncol=3)
	colnames(pathways_biased_true) <- c("traitfirst", "simultaneous", "climatefirst")
	# count truenumber of first transitions into % pathways
	for(i in 1:length(sims_traitfirst)) {
		cat("\nworking on tree", i)
		x <- sims_traitfirst[[i]]

		result <- countFirstTransitions(x)
		pathways <- (result[1:3]/result[4]) * 100  # [1:3] are the number of times a pathway is taken;  [4] is total number of first transitions
		pathways_biased_true <- rbind(pathways_biased_true, pathways)
	}

	# ----------------------
	# --  Checkpoint 11
	# ----------------------
	save(file="./intermediate_results/11_results_simulations.Rdata", list=c("pathways_biased_true", "pathways_biased_inf", "pathways_true", "pathways_inf"))

	#########################
	#  task 4.  plot results.
	#########################
	rm(list=ls())

	load("./intermediate_results/11_results_simulations.Rdata")

	pdf(width=8.9/2.54, height=8.9/2.54, file="./figures/Figure_simulations.pdf")

	par(mar=c(1.5,1.5,0.5, 0.5))
	simulcex <- 0.5
	plot(NULL, bty="o",
		xlim=c(0,100), ylim=c(0,100), 
		type="n", xaxt="n", yaxt="n")
		axis(side=1, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=-3, tck=-0.02)
		mtext("True % trait-first pathway (counted from simulated history)", cex=0.6, side=1, line=0.5)

	axis(side=2, at=seq(0, 100, 25), labels=paste(seq(0, 100, 25), "", sep=""),
		lwd.ticks=0.5, cex.axis=simulcex, padj=2, tck=-0.02)
	mtext("Inferred % trait-first pathway", cex=0.6, side=2, line=0.8)
	x <- pathways_biased_true[,1]
	y <- pathways_biased_inf[,1]
	col <- rep("grey20", 200)
	col[pathways_biased_inf[,1]<56] <- "red"
	points(x, y, col=col, pch=19, cex=0.666)
	x <- pathways_true[,1]
	y <- pathways_inf[,1]
	col <- rep("grey20", 200)
	col[pathways_inf[,1]>56 | pathways_inf[,1]<44] <-"red"
	points(x, y, col=col, pch=19, cex=0.666)

	x <- 0
	y <- 95
	points(x, y, col="grey20", pch=19, cex=1)
	text(x, y=y+2, labels="True history inferred", pos=4, cex=0.6)
	text(x, y=y-2, labels="as \"far more\" likely", pos=4, cex=0.6)
	y <- y-10
	points(x, y, col="red", pch=19, cex=1)
	text(x, y=y+2, labels="Incorrect history inferred", pos=4, cex=0.6)
	text(x, y=y-2, labels="as \"far more\" likely", pos=4, cex=0.6)

	mtext("B", side=1, line=0, adj=-2/30)
	
	dev.off()

	# conclusion: see Figure_simulations.pdf:  the method performs qualitatively fairly well when there is an underlying, much more likely trajectory (right column of dots). However, when both trajectories are equally likely, the method typically (in 154 out of 200 simulated histories) infers a higher contribution for a trajectory nevertheless (i.e., false positives, or 77% type I error).



