
Description of files (representing data, scripts, and results) associated with a critique of Zanne AE et al. (2014) Three keys to the radiation of angiosperms into freezing environments. Nature 506: 89–92.


################################

The folder in which the README is located contains the following files:
 - README.md  			 :	This very file		
 - annotated_analyses_script.R  :  The file with all annotated R scripts to repeat our 
 							analyses. It demonstrates the implementation of our analyses,
 							but running them serially takes very long.  Therefore, we 
 							actually executed them on the OSCAR computer cluster at the 
 							Centre for Computation and Visualization at Brown University. 
 							To facilitate following all our steps, we provide saves of 
 							intermediate results after each section of the script.
 - R_source_custom_functions.R  :  This file contains custom R functions which we wrote to 
 							facilitate the analyses.  It is read from the R script.
 - results_ClimateScoringSensitivity.txt  :  This is a tab-delimited file containing the 
 							results of analyses that differed in the handling of the 
 							climate data (annotated_analyses_script.R contains the
 							relevant annotations). It contains the following elements 
 							indicated by headers:
		file, treetype	treename : descriptors of the analysis
		ntaxa	: number of taxa in the tree
		n00		: number of taxa freezing unexposed without freezing-adaptated trait
		n01		: number of taxa freezing unexposed with freezing-adaptated trait
		n10		: number of taxa freezing exposed without freezing-adaptated trait
		n11		: number of taxa freezing exposed with freezing-adaptated trait
		groom	: imprecise and doubtful coordinates from GbIF data excluded?
		minclimdat : minimum number of coordinates required for including a species 
		perc	: percentile threshold to score a species as freezingExposed
		trait	: vessel, phenology, or woodiness data set	
		ratemat	: were simultaneous transitions ok (simultOK) or not (NULL)
		loglik	: log likelihood
		AICc	: AICc score
		trans_00.01	  : transition rate from state 00 to state 01
			...	
		trans_11.10	  : transition rate from state 11 to state 10	
		traitfirst	  : % contribution of trait-first pathway
		simult		  : % contribution of simultaneous pathway
		climfirst	  : % contribution of climate-first pathway
		climratesum	  : sum of all climate-state transition rates
		traitratesum  : sum of all trait-state transition rates	
		clim_vs_trait_labil	: 
		persistence00 : persistence times state 00	
			...	
		persistence11 : persistence times state 11
		diagnosis	  : did the analyses execute correctly?
		analysis	  : analysis identifier
		dAIC		  : delta-AIC score for preference of the model used
		

And it contains the following subfolders:
 - ./dryad/  			 :  Contains the data as provided by Zanne et al.
 - ./handling_climate_data/  :  Contains data, scripts, and results, pertaining to 
 							summarizing gbif localities into a single minimum-temperature
 							value per species
 - ./intermediate_results/  :  Contains intermediate results, that is, results generated
 							and processed in R at 11 intermediate stages of the analyses.
 - ./figures/				:  Contains additional figures generated using the R script.
 							annotated_analyses_script.R
 

################################

The subfolder "./dryad/" contains 5 files, that are copied unaltered from the datadryad digital repository where Zanne et al. uploaded their data
 - GlobalLeafPhenologyDatabase.csv		: leaf phenology data (abbreviated "phenol")
 - GlobalVessel_A_data.txt				: conduit size data (abbreviated "vessel")
 - GlobalWoodinessDatabase.csv			: woodiness data (abbreviated "woody")
 - MinimumFreezingExposure.csv			: climate data as handled by Zanne et al.
 - Vascular_Plants_rooted.dated.tre		: Zanne et al.'s phylogeny


################################

The subfolder "./handling_climate_data/" contains
 - getTminPerSpecies.py  :  An annotated phython script to pars the gbif locality data
 							as we received it from Zanne et al. and summarize it as one
 							minimal value per species.
 - allHerbaria_ADM1_badCoords.txt  :  The coordinates that represent the corners of the 
 							0.01° x 0.01° gridcells that contain centroids of major 
 							political areas and herbaria that house > 1,500,000 specimens. 
 							See annotation in phython script for details.
 - gbif_cords_tmin_species-summ.txt  :  The extensive output file of the python script.  
 							It contains more information that just for the R analyses.
 - errors.out			 :  The error outfile containing the lines in the gbif locality 
 							data file that needed to be excluded due to presence of NA.
 - climate_data_grooming_levels_JdV.txt  :  The slimmed version of the file 
 							gbif_cords_tmin_species-summ.txt, containing only those 
 							columns that contain data analyzed in the downstream analyses.
 							NB: the file with the actual gbif coordinate data is not
 							provided because it is 2.25 GB in size. We received from Zanne
 							et al.
  
 
################################
 
The subfolder "./intermediate_results/" contains intermediate results, listed and named in order of appearance in the R script. 
 
  
################################
 
The subfolder "./figures/" contains two figures generated by the R script, one comparable to Figure 1 in our critique (provided to facilitate comparison between the analysis settings and the individual results), the other illustrates the results of the simulation study.
  
  
  
  
  
  