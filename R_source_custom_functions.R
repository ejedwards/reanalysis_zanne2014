### R functions used in annotated_analyses_script.R



	# makeTransparent() is a function to easily add alpha to a color
	# from http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
	makeTransparent <- function(someColor, alpha=50){  # pass alpha on 0-255 scale
		newColor<-col2rgb(someColor)
  		apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    	blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
	} 
	  
	# runner() is a function to transforms simmap tree object x into a tree and data
	# object that can be analyzed with corDISC
	# The function assumes that along the tree four characters have been simulated 
	# representing states 00, 01, 10 and 10 
	runner <- function(x, rate_mat=FALSE) { 
		if(rate_mat){
			# make an index matrix to allow simultaneous transitions	
			rate_mat <- matrix(data=c(NA, 1, 2, 9,  3, NA, 10, 4,  5, 11, NA, 6,  12, 7, 8, NA), byrow=FALSE, nrow=4)
			rownames(rate_mat) <- c("(0,0)", "(0,1)", "(1,0)", "(1,1)")
			colnames(rate_mat) <- rownames(rate_mat)
		} else {rate_mat=NULL}
		# get data in shape
			# there are  four simulated states, 1:4 
			# these are interpreted as 00 (ancestral), 01 (only trait shifted), 10 (only climate shifted), 11 (derived), respectively.
			simdat <- data.frame(spec=x$tip.label, clim=NA, trait=NA)
			simdat[,"clim"][x$states%in%c(1,2)] <- "0"
			simdat[,"clim"][x$states%in%c(3,4)] <- "1"
			simdat[,"trait"][x$states%in%c(1,3)] <- "0"
			simdat[,"trait"][x$states%in%c(2,4)] <- "1"
		# do the analysis
			simres <- corDISC(x, simdat, rate.mat=rate_mat, model="ARD")
		# return the object
			return(simres)
	}

	# getPathways() is a function that takes as input a 4x4 transition rate matrix and 
	# calculates from it the % contribution of the different pathways from (0,0) to (1,1),
	# as clarified to us upon contact with Jeremy Beaulieu and Brian O'Meara around 
	# May 19, 2014
	getPathways <- function(pars) { 
		Q <- matrix(,6,6)
		Q[1:3, 1:3] <- pars[1:3, 1:3]
		Q[2,4] <- pars[2,4]  # 01-11 i.e. first trait then climate
		Q[1,5] <- pars[1,4]  # 00-11 i.e. simultaneous
		Q[3,6] <- pars[3,4]  # 10-11 i.e. first climate then trait
		Q[is.na(Q)] <- 0
		diag(Q) <- 0
		diag(Q) <- -rowSums(Q)
		liks <- c(1,0,0,0,0,0)
		pathways <- (expm::expm(t(Q) * 140*10^6, method=c("Ward77")) %*% liks)[4:6]  # final cols only
		pathways <- matrix(pathways*100, nrow=1)
		colnames(pathways) <- c("traitfirst", "simultaneous", "climfirst")

		pathways
	}


	# countFirstTransitions() is a function to count from a simmap trees the number of 
	# times the state (1,1) (i.e. final state) was reached from states (0,0), (0,1), and
	# (1,0). It only counts the first transitions (looking from the root forward in 
	# time), not the total number of transitions to (1,1).
	countFirstTransitions <- function(simmapTree, final="4") {
		x <- simmapTree
		# first round: paint downstream from first transitions to final state to obscure nested transitions
		for(i in 1:length(x$maps)) {
			state4onbranch <- names(x$maps[[i]]) %in% final   # length of element gives # transitions, sum > 1 if "4" present
			if(length(state4onbranch)>1){  # i is a transition branch
				if(sum(state4onbranch[2:length(state4onbranch)])>0) {# i contains transition to 0
					# paint nodes downstream into state "4"
					# get descendent of branch i, x$edge[i,2] to x
					desc <- x$edge[i,2]
					if (desc > length(x$tip.label)) { # its not a branch to a tip
						# paint clade from desc down to state "4"
						x <- paintSubTree(x, node=desc, state=final, stem=FALSE)
					} else {} # do nothing
				}
			} else {} # i is not a transition branch
		}	
		# second round: count branches with transitions to final state
		from1 <- 0  # state 1 can be interpreted as 00; from 1 is simultaneous
		from2 <- 0  # state 2 can be interpreted as 01; from 2 is traitfirst
		from3 <- 0  # state 3 can be interpreted as 10; from 3 is climatefirst
		for(i in 1:length(x$maps)) {
			state4onbranch <- names(x$maps[[i]]) %in% final
		
			if( length(state4onbranch)>1 & sum(state4onbranch[2:length(state4onbranch)])>0 ){
				# count if its a transition branch & its a transition to final state (default "4")
				# infer what the pre-final state was
				fromwhere <- names( x$maps[[i]] )[ which(state4onbranch==TRUE)[1] - 1 ]
				if(fromwhere=="1"){from1 <- from1 + 1}
				if(fromwhere=="2"){from2 <- from2 + 1}
				if(fromwhere=="3"){from3 <- from3 + 1}
			}
		}
		# final step: join result
		result <- c(from2, from1, from3, sum(from1, from2, from3))
		names(result) <- c("traitfirst", "simult", "climfirst", "total_first")
		return(result)
	}


	# fileListToResultTable() is a function to take a list of files that contain corDISC 
	# results, load them one by one, and summarize them in a table which is written to a 
	# file
	fileListToResultTable <- function(fileList, outputfile="../output.txt", append=FALSE) {
	cat("\nNumber of files:", length(fileList), "\n")	
	
	if(!(append)) { # start new file
		cat("file\ttreetype\ttreename\tntaxa\tn00\tn01\tn10\tn11\tgroom\tminclimdat\tperc\ttrait\tratemat\tloglik\tAICc\ttrans_00-01\ttrans_00-10\ttrans_00-11\ttrans_01-00\ttrans_01-10\ttrans_01-11\ttrans_10-00\ttrans_10-01\ttrans_10-11\ttrans_11-00\ttrans_11-01\ttrans_11-10\ttraitfirst\tsimult\tclimfirst\tclimratesum\ttraitratesum\tclim_vs_trait_labil\tpersistence00\tpersistence01\tpersistence10\tpersistence11\tdiagnosis", file=outputfile, append=append)
	} else {cat("\n", file=outputfile, append=append)}
	for(i in 1:length(fileList)) {
		cat("\n", i)
		inputfile <- fileList[[i]]
	
		# clean up previous analyses
		if("result" %in% ls()) {rm("result")}

	### fix filename
		filename <- sub("rate.mat_all", "simultOK", inputfile, fixed=TRUE)


	### extract elements from file name
		elem <- strsplit(filename, split="_")[[1]]  # split name into its elements
		x <- length(elem)
		ratemat <- elem[x-6]
		treename <- elem[x-5]
		if (elem[x-5] %in% c("asterids", "rosids", "magnolids", "resttree", "monocots")) {
			treetype <- "compound"
		} else { treetype <- elem[x-5] }
		groom <- elem[x-3]
		minclimdat <- elem[x-2]
		perc <- elem[x-1]
		trait <- strsplit(elem[x], ".", fixed=TRUE)[[1]][1]  # remove Rdata 
		
	### extract elements from the result
		load(inputfile) # gives element result (or failed with a reason)
		
		if(exists("result")) {
			doStats <- TRUE
			diagnosis <- "ok"
		} else {
			cat(paste("file", inputfile ,"does not contain object called result"))
			diagnosis <- "notOk" 
		}	
	
		if(doStats){
			ntaxa <- length(result$phy$tip.label)
			loglik <- result$loglik
			AICc <- result$AICc
			pars <- result$solution
			data <- result$data
		
			# get pathways
			pathways <- getPathways(pars)
			traitfirst <- pathways[1]
			simult <- pathways[2]
			climfirst <- pathways[3]
			
			# lability ratio:  
			   # sum clim transitions / sum trait transitions
			   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
			   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
			   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
			   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
			climratesum  <- sum(pars[1,3], pars[2,4], pars[1,4], pars[3,1], pars[4,2], pars[4,1], na.rm=TRUE)
			traitratesum <- sum(pars[1,2], pars[3,4], pars[1,4], pars[2,1], pars[4,3], pars[4,1], na.rm=TRUE)
			clim_vs_trait_labil <- climratesum / traitratesum
		
			# character data
			char00 <- sum(data[,1]=="0" & data[,2]=="0")
			char01 <- sum(data[,1]=="0" & data[,2]=="1")
			char10 <- sum(data[,1]=="1" & data[,2]=="0")
			char11 <- sum(data[,1]=="1" & data[,2]=="1")
			
			# persistence times
			   # sum clim transitions / sum trait transitions
			   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
			   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
			   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
			   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
			persis00  <- 1/rowSums(pars, na.rm=TRUE)[1]
			persis01  <- 1/rowSums(pars, na.rm=TRUE)[2]
			persis10  <- 1/rowSums(pars, na.rm=TRUE)[3]
			persis11  <- 1/rowSums(pars, na.rm=TRUE)[4]
		} 
		if(!doStats) {
		
		### list of NAs
			ntaxa <- NA
			loglik <- NA
			AICc <- NA
			pars <- matrix(NA,4,4)
			traitfirst <- NA
			simult <- NA
			climfirst <- NA
			climratesum <- NA
			traitratesum <- NA
			clim_vs_trait_labil <- NA
			char00 <- NA
			char01 <- NA
			char10 <- NA
			char11 <- NA
			persis00  <- NA
			persis01  <- NA
			persis10  <- NA
			persis11  <- NA
		}
	
	### write all elements to a line in the file
		
		cat("\n", file=outputfile, append=TRUE)
		
		cat(filename, treetype, treename, ntaxa, char00, char01, char10, char11, groom, minclimdat, perc, trait, ratemat, loglik, AICc, pars[1,2], pars[1,3], pars[1,4], pars[2,1], pars[2,3], pars[2,4], pars[3,1], pars[3,2], pars[3,4], pars[4,1], pars[4,2], pars[4,3], traitfirst, simult, climfirst, climratesum, traitratesum, clim_vs_trait_labil, persis00, persis01, persis10, persis11, diagnosis, sep="\t", file=outputfile, append=TRUE)
	}
}


	# function resultTableToSelectedModels() reads an object created by function 
	# fileListToResultTable() and returns similar table file but with only best models 
	# for each analysis type. 
	resultTableToSelectedModels <- function( 
		res, 
		treetype="tree", testtype="perSub", 
		subset=FALSE, trait=NULL, 
		outputfile=NULL, append=FALSE
	)
	# res is the input generated by fileListToResultTable()
	# if treetype="compound" : get compound weighted-average transition rates
	# if testtype="jointly"  : do model selection per group of 5 
	#						   (all with or all without simultaneous transitions)
	# if testtype="direct"   : simply take the weighted average of the trajectories, 
	#						   which is the way of the paper.
	# if subset=TRUE		 : use trait and treetype to select relevant subset of results
	# if outputfile specified, write results to a file
	{
		# subset input; defaults to FALSE leading to taking everything
		if(subset == TRUE) {
			res <- res[res$trait==trait & res$treetype==treetype, ]
		}

		# compound results when multiple analyses were run across a tree as for the woody data
		# if so, select option to calculate pathways. possibilities: perSub, jointly, direct
		# [NB we are assuming each partition gets its own parameters, as in Zanne et al.]

		if(treetype=="compound" & testtype=="direct") {   # simple numeric summary of the rates

			res$analysis <- with(res, paste(trait, groom, minclimdat, perc))

			# do joint model selection 

			res$analysis <- with(res, paste(trait, groom, minclimdat, perc))
			res$dAIC <- NA			
			res2 <- matrix(NA, nrow=0, ncol=dim(res)[2])
			colnames(res2) <- colnames(res)

			for (thisAnalysis in unique(res$analysis)) {	# for each of the "super"analyses of 5 subtrees...
				# grab analyses belonging together
				rows_to_use <- which(res$analysis == thisAnalysis)
				res_tmp <- res[ rows_to_use, ]
				
				# grab likelihood product under ratemat and NULL and compare AIC				
				lik_NULL <-	sum(res_tmp$loglik[res_tmp$ratemat == "NULL"])
				lik_ratemat <-	sum(res_tmp$loglik[res_tmp$ratemat == "ratematall"])
				
				AIC_NULL    <- 2*5*8 - 2*lik_NULL
				AIC_ratemat <- 2*5*12 - 2*lik_ratemat
				dAIC <- abs(AIC_NULL-AIC_ratemat)
				if(AIC_NULL < AIC_ratemat) {best <- "NULL"} else {best <- "ratematall"}
				
				# grab the relevant subtrees for the weighted-average transition rates
				res_tmp <- res_tmp[ res_tmp$ratemat==best & res_tmp$treename != "resttree", ]

				pathways_all <- res_tmp[, c("traitfirst", "simult", "climfirst")]
				specs <- res_tmp$ntaxa
				traitfirst <- sum(pathways_all[,1]*(specs/sum(specs)))
				simult <- sum(pathways_all[,2]*(specs/sum(specs)))
				climfirst <- sum(pathways_all[,3]*(specs/sum(specs)))

				# we pick the shared elements from the these and subsititue the rest 
				x <- res_tmp[1, ]
				outputline <- c(x[[1]], x[[2]], NA, x[[4]], x[[5]], x[[6]], x[[7]], x[[8]], x[[9]], x[[10]], x[[11]], x[[12]], NA, NA, NA, rep(NA, 12), traitfirst, simult, climfirst, rep(NA, 7), "compound", thisAnalysis, dAIC)
				outputline <- matrix(data=outputline, nrow=1)
				res2 <- rbind(res2, outputline)  # add it to the result matrix
			}
			res2 <- data.frame(res2, stringsAsFactors=FALSE)	
		}
		if(treetype=="compound" & testtype=="perSub") {  

			res$analysis <- with(res, paste(trait, groom, minclimdat, perc))

			# do modelselection per subtree
			# (perhaps this is not how it's done and need to grab all simult or all non-simult models)
			analysis_tmp <- with(res, paste(trait, groom, minclimdat, perc, treename))
			res_tmp <- matrix(NA, nrow=0, ncol=dim(res)[2]+1)
			colnames(res_tmp) <- c(colnames(res), "dAIC")

			for(i in unique(analysis_tmp)) { # do model selection per subtree
				tmp <- res[analysis_tmp==i,] # grab models for an analysis
				best <- tmp[which.min(as.numeric(tmp$AIC)), ] # select obtimal one 
				best$dAIC <- abs(as.numeric(tmp$AIC[1])-as.numeric(tmp$AIC[2])) # get dAIC for that one
				res_tmp <- rbind(res_tmp, best)  # add it to the result matrix
			}
			
			res2 <- matrix(NA, nrow=0, ncol=dim(res_tmp)[2])
			colnames(res2) <- colnames(res_tmp)

			#  get weighted average results using best subtrees  
			for (thisAnalysis in unique(res$analysis)) {	# for each of the "super"analyses of 5 subtrees...

				# grab analyses belonging together
				rows_to_use <- which(res_tmp$analysis == thisAnalysis)
				cols_to_use <- which(grepl("trans_", colnames(res_tmp)))

				# get the weighted-average transition rates
				res_tmp[rows_to_use,cols_to_use][is.na(res_tmp[rows_to_use,cols_to_use])] <- 0
				tmp <- colSums(res_tmp[rows_to_use,cols_to_use] * res_tmp[rows_to_use,"ntaxa"])/	sum(res_tmp[rows_to_use,"ntaxa"])
				weightTransMat <- matrix(c(NA, tmp[1:4], NA, tmp[5:8], NA, tmp[9:12], NA), nrow=4, byrow=TRUE)
				pars <- weightTransMat

				# get pathways
				pathways <- getPathways(pars)
				traitfirst <- pathways[1]
				simult <- pathways[2]
				climfirst <- pathways[3]
			
				# get lability ratio:  
				   # sum clim transitions / sum trait transitions
				   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
				   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
				   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
				   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
				climratesum  <- sum(pars[1,3], pars[2,4], pars[1,4], pars[3,1], pars[4,2], pars[4,1], na.rm=TRUE)
				traitratesum <- sum(pars[1,2], pars[3,4], pars[1,4], pars[2,1], pars[4,3], pars[4,1], na.rm=TRUE)
				clim_vs_trait_labil <- climratesum / traitratesum
			
				# get persistence times
				persis00  <- 1/rowSums(pars, na.rm=TRUE)[1]
				persis01  <- 1/rowSums(pars, na.rm=TRUE)[2]
				persis10  <- 1/rowSums(pars, na.rm=TRUE)[3]
				persis11  <- 1/rowSums(pars, na.rm=TRUE)[4]
	
				# write line to the output table
				# corresponding elements in res_tmp were: 
				#		cat(filename, treetype, treename, ntaxa, char00, char01, char10, char11, groom, minclimdat, perc, trait, ratemat, loglik, AICc, pars[1,2], pars[1,3], pars[1,4], pars[2,1], pars[2,3], pars[2,4], pars[3,1], pars[3,2], pars[3,4], pars[4,1], pars[4,2], pars[4,3], traitfirst, simult, climfirst, climratesum, traitratesum, clim_vs_trait_labil, persis00, persis01, persis10, persis11, diagnosis, sep="\t", file=outputfile, append=TRUE)
				# we pick the shared elements from the these and subsititue the rest 
				x <- res_tmp[rows_to_use[1], ]
				outputline <- c(x[[1]], x[[2]], NA, x[[4]], x[[5]], x[[6]], x[[7]], x[[8]], x[[9]], x[[10]], x[[11]], x[[12]], NA, NA, NA, pars[1,2], pars[1,3], pars[1,4], pars[2,1], pars[2,3], pars[2,4], pars[3,1], pars[3,2], pars[3,4], pars[4,1], pars[4,2], pars[4,3], traitfirst, simult, climfirst, climratesum, traitratesum, clim_vs_trait_labil, persis00, persis01, persis10, persis11, "compound", thisAnalysis, NA)
				outputline <- matrix(data=outputline, nrow=1)
				res2 <- rbind(res2, outputline)  # add it to the result matrix
			}
			res2 <- data.frame(res2, stringsAsFactors=FALSE)	
		}
		
		if(treetype=="compound" & testtype=="jointly") {  
		# do model selection for the 5 subtrees jointly
			
			res$analysis <- with(res, paste(trait, groom, minclimdat, perc))
			res$dAIC <- NA			
			res2 <- matrix(NA, nrow=0, ncol=dim(res)[2])
			colnames(res2) <- colnames(res)

			for (thisAnalysis in unique(res$analysis)) {	# for each of the "super"analyses of 5 subtrees...
#thisAnalysis <- unique(res$analysis)[1]
				# grab analyses belonging together
				rows_to_use <- which(res$analysis == thisAnalysis)
				res_tmp <- res[ rows_to_use, ]
				
				# grab likelihood product under ratemat and NULL and compare AIC				
				lik_NULL <-	sum(res_tmp$loglik[res_tmp$ratemat == "NULL"])
				lik_ratemat <-	sum(res_tmp$loglik[res_tmp$ratemat == "ratematall"])
				
				AIC_NULL    <- 2*5*8 - 2*lik_NULL
				AIC_ratemat <- 2*5*12 - 2*lik_ratemat
				dAIC <- abs(AIC_NULL-AIC_ratemat)
				if(AIC_NULL < AIC_ratemat) {best <- "NULL"} else {best <- "ratematall"}
				
				# grab the relevant subtrees for the weighted-average transition rates
#save.image("~/Desktop/tmp_dump.Rdata")				
				res_tmp <- res_tmp[ res_tmp$ratemat==best, ]
				cols_to_use <- which(grepl("trans_", colnames(res_tmp)))
				
				res_tmp[,cols_to_use][is.na(res_tmp[,cols_to_use])] <- 0 # because colSums doesn't like  NA

				tmp <- colSums(res_tmp[ ,cols_to_use] * res_tmp[,"ntaxa"])/	sum(res_tmp[,"ntaxa"])
				weightTransMat <- matrix(c(NA, tmp[1:4], NA, tmp[5:8], NA, tmp[9:12], NA), nrow=4, byrow=TRUE)
				pars <- weightTransMat


				# get pathways
				pathways <- getPathways(pars)
				traitfirst <- pathways[1]
				simult <- pathways[2]
				climfirst <- pathways[3]
			
				# get lability ratio:  
				   # sum clim transitions / sum trait transitions
				   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
				   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
				   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
				   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
				climratesum  <- sum(pars[1,3], pars[2,4], pars[1,4], pars[3,1], pars[4,2], pars[4,1], na.rm=TRUE)
				traitratesum <- sum(pars[1,2], pars[3,4], pars[1,4], pars[2,1], pars[4,3], pars[4,1], na.rm=TRUE)
				clim_vs_trait_labil <- climratesum / traitratesum
			
				# get persistence times
				persis00  <- 1/rowSums(pars, na.rm=TRUE)[1]
				persis01  <- 1/rowSums(pars, na.rm=TRUE)[2]
				persis10  <- 1/rowSums(pars, na.rm=TRUE)[3]
				persis11  <- 1/rowSums(pars, na.rm=TRUE)[4]
	
				# write line to the output table
				# corresponding elements in res_tmp were: 
				#		cat(filename, treetype, treename, ntaxa, char00, char01, char10, char11, groom, minclimdat, perc, trait, ratemat, loglik, AICc, pars[1,2], pars[1,3], pars[1,4], pars[2,1], pars[2,3], pars[2,4], pars[3,1], pars[3,2], pars[3,4], pars[4,1], pars[4,2], pars[4,3], traitfirst, simult, climfirst, climratesum, traitratesum, clim_vs_trait_labil, persis00, persis01, persis10, persis11, diagnosis, sep="\t", file=outputfile, append=TRUE)
				# we pick the shared elements from the these and subsititue the rest 
				x <- res_tmp[1, ]
				outputline <- c(x[[1]], x[[2]], NA, x[[4]], x[[5]], x[[6]], x[[7]], x[[8]], x[[9]], x[[10]], x[[11]], x[[12]], NA, NA, NA, pars[1,2], pars[1,3], pars[1,4], pars[2,1], pars[2,3], pars[2,4], pars[3,1], pars[3,2], pars[3,4], pars[4,1], pars[4,2], pars[4,3], traitfirst, simult, climfirst, climratesum, traitratesum, clim_vs_trait_labil, persis00, persis01, persis10, persis11, "compound", thisAnalysis, dAIC)
				outputline <- matrix(data=outputline, nrow=1)
				res2 <- rbind(res2, outputline)  # add it to the result matrix
			}
			res2 <- data.frame(res2, stringsAsFactors=FALSE)	
		}
	
		
		if(treetype != "compound") {	
			# do modelselection; keep only best per analysis
			res$analysis <- with(res, paste(trait, groom, minclimdat, perc))
			res2 <- matrix(NA, nrow=0, ncol=dim(res)[2]+1)
			colnames(res2) <- c(colnames(res), "dAIC")
			for(i in unique(res$analysis)) { # do model selection
				tmp <- res[res$analysis==i,] # grab models for an analysis
				best <- tmp[which.min(as.numeric(tmp$AIC)), ] # select obtimal one 
				best$dAIC <- abs(as.numeric(tmp$AIC[1])-as.numeric(tmp$AIC[2])) # get dAIC for that one
				res2 <- rbind(res2, best)  # add it to the result matrix
			}
			res2 <- data.frame(res2, stringsAsFactors=FALSE)	
		}
		
		# return result and write it down 
		if(!is.null(outputfile)) {
			if(append==FALSE) {  # write the header; klopt
				header <- "file\ttreetype\ttreename\tntaxa\tn00\tn01\tn10\tn11\tgroom\tminclimdat\tperc\ttrait\tratemat\tloglik\tAICc\ttrans_00.01\ttrans_00.10\ttrans_00.11\ttrans_01.00\ttrans_01.10\ttrans_01.11\ttrans_10.00\ttrans_10.01\ttrans_10.11\ttrans_11.00\ttrans_11.01\ttrans_11.10\ttraitfirst\tsimult\tclimfirst\tclimratesum\ttraitratesum\tclim_vs_trait_labil\tpersistence00\tpersistence01\tpersistence10\tpersistence11\tdiagnosis\tanalysis\tdAIC\n"
			cat(header, file=outputfile, append=FALSE)
			}
			for (i in 1:dim(res2)[1]) {
				cat(paste(res2[i, ]), sep="\t", file=outputfile, append=TRUE)
				cat("\n",   file=outputfile, append=TRUE)
			}
		}
		return(res2)
	}


	# pickSimul() is a function to simply compare two corDISC result objects based on $AIC 
	# element to print which one is better fitting
	pickSimul <- function(nos, s) {
		min    <- which.min( c(nos$AIC, s$AIC) )
		result <- paste("Simultaneous changes", c("not", "")[min], "supported; use model", c("null", "ratemat")[min])
		cat(result)
	}

	# getPathways_print() is a function that takes as input a corDISC result object and 
	# uses it to print the % pathway and other results (such as trait lability and 
	# persistence times)
	getPathways_print <- function(corDISCobject){
		result <- corDISCobject
			
		ntaxa <- length(result$phy$tip.label)
		pars <- result$solution
		data <- result$data
	
		# get pathways
		pathways <- getPathways(pars)
		traitfirst <- pathways[1]
		simult <- pathways[2]
		climfirst <- pathways[3]
		
		# lability ratio:  
		   # sum clim transitions / sum trait transitions
		   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
		   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
		   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
		   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
		climratesum  <- sum(pars[1,3], pars[2,4], pars[1,4], pars[3,1], pars[4,2], pars[4,1], na.rm=TRUE)
		traitratesum <- sum(pars[1,2], pars[3,4], pars[1,4], pars[2,1], pars[4,3], pars[4,1], na.rm=TRUE)
		clim_vs_trait_labil <- climratesum / traitratesum
	
		# character data
		char00 <- sum(data[,1]=="0" & data[,2]=="0")
		char01 <- sum(data[,1]=="0" & data[,2]=="1")
		char10 <- sum(data[,1]=="1" & data[,2]=="0")
		char11 <- sum(data[,1]=="1" & data[,2]=="1")
		
		# persistence times
		   # sum clim transitions / sum trait transitions
		   #  clim  transitions are 00-10, 01-11, 00-11; reverse 10-00, 11-01, 11-00
		   #			in Q:		[1,3], [2,4], [1,4];         [3,1], [4,2], [4,1]
		   #  trait transitions are 00-01, 10-11, 00-11; reverse 01-00, 11-10, 11-00
		   #			in Q:		[1,2], [3,4], [1,4];         [2,1], [4,3], [4,1]
		persis00  <- 1/rowSums(pars, na.rm=TRUE)[1]
		persis01  <- 1/rowSums(pars, na.rm=TRUE)[2]
		persis10  <- 1/rowSums(pars, na.rm=TRUE)[3]
		persis11  <- 1/rowSums(pars, na.rm=TRUE)[4]
	
	
		phrase <- paste(
			"Summary of corDISC result:",
			"\n--------------------------",
			"\npercent trait-first  : ", round(traitfirst, 2),
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
			"\nnumber taxa with clim.ances & trait.ances: ", char00,
			"\nnumber taxa with clim.ances & trait.deriv: ", char01,
			"\nnumber taxa with clim.deriv & trait.ances: ", char10,
			"\nnumber taxa with clim.deriv & trait.deriv: ", char11,
			"\n",
			"\nfraction of taxa with freezing-adapted trait state: ", round((char01+char11)/ntaxa, 3),
			sep=""
		)
		cat(phrase)
	} 





