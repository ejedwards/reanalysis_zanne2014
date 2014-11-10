#!/usr/bin/env python
import re
import time
import scipy.stats
# define function mean and median
def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
    return sorts[length / 2]


Debugging = False


# *************************************   outline   ************************************* #
# This script reads an input file (as we received it from Zanne et al.) containing 
# Gbif records of latitude, longitude, and the associated minimum annual temperature 
# for known species. 
# 
# Format of the input file:
# "Scientific_name_interpreted","Longitude","Latitude","tmin"
# "Abarema adenophora",-77.1333,3.93333,22
# "Abarema adenophora",-72.1333333,0.65,21
# "Abarema adenophora",-71.1625833,2.1441111,21
# "Abarema adenophora",-68.0588889,4.5077778,22.5
# "Abarema adenophora",-68.0588889,4.5077778,22.5
#
# The objective of the script is to summarize these records to a single value per species.
# Zanne et al. did this by obtaining the value of Tmin at the 2.5 percentile, to subsequently check whether it is larger than 0°C, or equal or smaller than 0°C.  We extend this approach in three ways. 
# (1) We impose a minimal sampling size (1, 3, or 10 records per species)
# (2) We groom the data more extensively, by eliminating duplicate coordinates (set 
# _noDup), or in addition eliminating implausible and imprecise records (set _groomed). 
# with implausible we mean those records exactly at the centroid of a country or primary 
# administrative division, e.g. a state or providence, and those records exactly at the 
# location of a herbarium that houses >1,500,000 specimens according to 
# http://en.wikipedia.org/wiki/List_of_herbaria (version of 18 February 2014). We define 
# "exactly" as occurring anywhere within the 0.01° x 0.01° gridcell of the herbarium or 
# political centroid. Imprecise records are defined as those records that have only a 
# single digit, because those are less precise than the climate data. 
# (3) Because many species have localities that experience freezing and those that do not, 
# we calculate the species' value of Tmin by using the score at alternative percentiles. 
# We used percentiles 0 (i.e. if a species has a single record ≤0°C, it is freezing 
# exposed), 2.5, 10, 50 (median), 90, 97.5, 100 (i.e. only when all records of a species 
# experience freezing, the species is freezing exposed).
#
# The script writes the extensive results in a file "gbif_cords_tmin_species-summ.txt", 
# which is provided.  The downstream analyses in R use a version of this data in 
# which we manually removed unnecessary columns using Excel.
#
# outline of the program:
# for each line of the csv value, except header)
#	split it, remember as tuple
# 	check whether it contains NA
# 		if so:  
# 			write the line down in the error file and forget about it
# 		if not: 
# 			read species name from LineElements
# 			check whether different from previous species
# 				if same as previous:
# 					append whole line to running list of tuple
# 				if different species than previous species:
# 	
# 					* DETERMINE _all AND _noDup RECORDS from tuple using set()
#					* DETERMINE _groomed by comparison to "bad" coordinates
#						- write down counts to file
# 					* CALCULATE STATS ON _all, _noDup, _groomed records
# 						- write down lat lon stats to file
# 					* RESET CURRENT LISTS AND COUNTERS
# 					* APPEND LINE TO NEW running tuple
# 		remember which line we just read
#		every 10000 lines, print a summary of what has been achieved
# 		next line

# to keep track of speed:
TimeStart = time.time()

# define file names   ####### insert correct path for your local system before executing
InFileName 			= '/foo/bar/github/handling_climate_data/gbif_cords_tmin.csv'  
# 	this is the file with the Gbif locality data which we received from Zanne et al. It is 
#	not provided because it is very large (2gb)
PolitCoordFileName	= '/foo/bar/github/handling_climate_data/allHerbaria_ADM1_badCoords.txt'
# 	this file contains the four corners of each 0.01 x 0.01 gridcell that contains a herbarum or political coordinate 
BrokenLinesFileName	= '/foo/bar/github/handling_climate_data/errors.out'
#  Error file for badly formatted records
SpecSummFileName 	= '/foo/bar/github/handling_climate_data/gbif_cords_tmin_species-summ.txt'
#  Output file

# provide header of species summary file SpecSummFileName with trailing \n
SpecSummFileHeader = "Scientific_name_interpreted\tRecords_all\tRecords_noDup\tRecords_groomed\tLat_min_all\tLat_min_groomed\tLat_min_noDup\tLon_min_all\tLon_min_groomed\tLon_min_noDup\tLat_med_all\tLat_med_groomed\tLat_med_noDup\tLon_med_all\tLon_med_groomed\tLon_med_noDup\tLat_max_all\tLat_max_groomed\tLat_max_noDup\tLon_max_all\tLon_max_groomed\tLon_max_noDup\tTmin_all_score0perc\tTmin_groomed_score0perc\tTmin_noDup_score0perc\tTmin_all_percmin\tTmin_all_perc0010\tTmin_all_perc0025\tTmin_all_perc0050\tTmin_all_perc0100\tTmin_all_percmed\tTmin_all_perc0900\tTmin_all_perc0950\tTmin_all_perc0975\tTmin_all_perc0990\tTmin_all_percmax\tTmin_groomed_percmin\tTmin_groomed_perc0010\tTmin_groomed_perc0025\tTmin_groomed_perc0050\tTmin_groomed_perc0100\tTmin_groomed_percmed\tTmin_groomed_perc0900\tTmin_groomed_perc0950\tTmin_groomed_perc0975\tTmin_groomed_perc0990\tTmin_groomed_percmax\tTmin_noDup_percmin\tTmin_noDup_perc0010\tTmin_noDup_perc0025\tTmin_noDup_perc0050\tTmin_noDup_perc0100\tTmin_noDup_percmed\tTmin_noDup_perc0900\tTmin_noDup_perc0950\tTmin_noDup_perc0975\tTmin_noDup_perc0990\tTmin_noDup_percmax\n"

# set desired level of grooming 
DesiredCoordPrecision		= 2		# to throw out records that have less than this number of digits.

# initiate trackers
LineNumber      	 		= 0		# keeps track of line number
PrintNumber			 		= 0		# keeps track of iterations to only print every x generation
SpeciesNumber  		  		= 0		# keeps track of number of species processed
BrokenLinesNumber	 		= 0		# keeps track of lines with 'NA'
BeyondStart 				= False # switches to TRUE after reading header
BadLocalitySpeciesCounter	= 0 	# keeps track of species that lose all records after grooming
BadSpeciesCounter			= 0 	# keeps track of species with only non-precise localities

# initiate objects to store lines in
SpeciesLineList 		= []	# List to store all lines (as tuples) for a species
SpeciesThisLine  		= '""'	# first species; keeps track of species processed in the previous line
SpeciesPreviousLine 	= '""'	# first species; keeps track of species processed in the previous line
BadCoordinates			= []	# List of tuples of the "bad" coordinates that match major herbaria
BadSpeciesTot			= []	# total number of records for each bad species
BadSpeciesBad			= []	# number of bad localities for each bad species
BadSpeciesName			= []	# names of species with bad localities


# read set of coordinates that represent herbaria and political centroids
PolitCoordFile	= open(PolitCoordFileName, 'r')		# csv input: meaningless coordinates that need to be excluded (herbaria + political areas ADM1 level)
for Line in PolitCoordFile:
	Line        = Line.strip('\n')
	LineElements= Line.split(',')
	BadTuple 	= (str(LineElements[0]), str(LineElements[1]))
	BadCoordinates.append(BadTuple)
BadCoordinates = set(BadCoordinates) 	# to have the object in useful format
PolitCoordFile.close()				

# open pipes to files
InFile       	= open(InFileName, 'r')				# csv input
SpecSummFile 	= open(SpecSummFileName, 'w')		# output: species summary stats
BrokenLinesFile	= open(BrokenLinesFileName, 'w')	# output: broken lines

# parse lines in the InFile:
# read all lines of a species, storing them as a list of tuples with this format:
#    SpeciesLineList = [ ('"gen spec"', (lat, long), clim), ('"gen spec"', (lat, long), clim), ... ]
# 1. filter lines with NA
# 2. filter lines with non-unique records using set()
# 3. filter lines with imprecise records (those with 1 dec degree only) 
# 4. filter lines with meaningless coordinates
# NB: There are 2384 species with records that match to a location of one of the 30 major 
# herbaria with 1,500,000+ specimens. Usually there are only few problematic localities 
# per species (1-10).

#	Read and start processing lines
for Line in InFile:
	LineNumber += 1

	# core: allowing for exception to handle header and first line of first species
	if BeyondStart:	
		# read line and split into its elements
		Line         = Line.strip('\n')
		LineElements = Line.split(',')

		# check for NAs in line. If not present, process the line
		if not('NA' in LineElements):
			# extract species from line
			SpeciesThisLine = LineElements[0]

			# seems like trailing 0's were dropped by Zanne et al., e.g. 10.10°N, 31.18°W 
			# is represented as 10.1, 31.18. Thus number of dec.digits is not identical
			# to actual coordinate precision.
			# but since precision of lat should be that of long, precision is 
			# max (dec.digits(lat), dec.digits(long))
			# therefore, maximally one in 100 imprecise records is incorrectly excluded.

			# Up to now, the coordinates are strings. So we can count the number of digits 
			# beyond a dot
		
			#  get coordinate precision			
			LatPrecision = 0
			LongPrecision = 0
			MinImprecision = ['166666', '333333', '666666', '833333']
			#   Decimal coordinate starts with:
			#	166666  = 10 min
			#	333333  = 20 min
			#	666666  = 40 min
			#	833333  = 50 min
			# precision for latitude:
			LatLength  = str(LineElements[2])
			LatLength  = LatLength.split('.')
			if (len(LatLength) == 2):					 
				# then it was splitted at the . and thus it has decimals
				if LatLength[1][0:6] in MinImprecision:  
					# then its exactly precise at imprecise minutes
					LatPrecision = int(1)
				else: 
					LatPrecision = int(len(LatLength[1]))
			else:
				LatPrecision = int(0)
			# precision for longitude:
			LongLength  = str(LineElements[1])
			LongLength  = LongLength.split('.')
			if (len(LongLength) == 2):					  
				# then it was splitted at the . and thus it has decimals
				if LongLength[1][0:6] in MinImprecision:  
					# then its exactly precise at imprecise minutes
					LongPrecision = int(1)
				else: 
					LongPrecision = int(len(LongLength[1]))
			else:
				LongPrecision = int(0)
							
			# prepare LineTuple; coordinates now will be floats
			tmp = [float(LineElements[x]) for x in range(1,4)]
			LineTuple = (SpeciesThisLine, (tmp[1], tmp[0]), tmp[2], LatPrecision, LongPrecision, )
				
			# if current species is the same as previous one, append to running list
			if SpeciesThisLine == SpeciesPreviousLine:
				SpeciesLineList.append(LineTuple)

			else: # read a new species? summarize the one we just completed reading
				SpeciesNumber += 1

				# * DETERMINE _all _noDup and _groomed records from tuple
				SpeciesRecords_all    			= SpeciesLineList   # for clarity
				SpeciesRecords_noDup 			= set(SpeciesLineList)
				SpeciesRecords_groomed 			= set([Records for Records in SpeciesRecords_noDup \
					if (Records[3] >= DesiredCoordPrecision) or (Records[4] >= DesiredCoordPrecision)])

				SpeciesRecords_groomed = SpeciesRecords_noDup - set([Tuple for Tuple in SpeciesRecords_noDup for Coord in BadCoordinates if (str(Coord[0]) in str(Tuple[1]) and str(Coord[1]) in str(Tuple[1]))])
				# i.e., subtracting bad coordinates from noDup set of tuples 												
				# if the coordinate is a bad coordinate	
	
				# get number of records. Use this to impose sampling size criterion
				Len_all = len(SpeciesRecords_all)
				Len_noDup = len(SpeciesRecords_noDup)
				Len_groomed = len(SpeciesRecords_groomed)
	
				if (Len_groomed > 0):  # do this to avoid calculating min(x) when x is empty
					# * CALCULATE STATS ON  _all, _noDup, _groomed
					# these are the necessary elements
					# for x in SpeciesLineList:
					#	print x			# full : ('"Abarema adenophora"', (3.93333, -77.1333), 22.0, 2, 3)
					#	print x[0]		# spec.: "Abarema adenophora"
					#	print x[1]		# coord: (6.35, -76.4333333)
					#	print x[1][0]	# lat  : 6.35
					#	print x[1][1]	# long : -76.4333333
					#	print x[2]  	# clim : 22.0
					#	print x[3]		# lat precision  : 4
					#	print x[4]		# long precision : 4
					
					# Extract latitude, longitude, and tmin from records
					Lats_all	 = [ x[1][0] for x in SpeciesRecords_all     ]
					Longs_all	 = [ x[1][1] for x in SpeciesRecords_all     ]
					Tmins_all	 = [ x[2] 	 for x in SpeciesRecords_all     ]
					Lats_noDup	 = [ x[1][0] for x in SpeciesRecords_noDup   ]
					Longs_noDup  = [ x[1][1] for x in SpeciesRecords_noDup   ]
					Tmins_noDup  = [ x[2] 	 for x in SpeciesRecords_noDup   ]
					Lats_groomed = [ x[1][0] for x in SpeciesRecords_groomed ]
					Longs_groomed= [ x[1][1] for x in SpeciesRecords_groomed ]
					Tmins_groomed= [ x[2] 	 for x in SpeciesRecords_groomed ]
					
					# Assign extracted lat lon values to dictionary 
					StatsDict = {
						'Lats_all'	  	: Lats_all, 
						'Lats_groomed'	: Lats_groomed,
						'Lats_noDup' 	: Lats_noDup,
						'Longs_all'	  	: Longs_all,
						'Longs_groomed'	: Longs_groomed,
						'Longs_noDup'	: Longs_noDup
					}
					SortedKeys = sorted(StatsDict.keys())
	
					# Do the stats
					Mins    = [min(   StatsDict.get(x)) for x in SortedKeys]
					Medians = [median(StatsDict.get(x)) for x in SortedKeys]
					Maxs    = [max(   StatsDict.get(x)) for x in SortedKeys]
					# (results ordered saved in the order of SortedKeys)
					
					# Assign extracted tmin values to dictionary 
					StatsDict = {
						'Tmins_all'	  	: Tmins_all,
						'Tmins_groomed'	: Tmins_groomed,
						'Tmins_noDup'	: Tmins_noDup
					}
					SortedKeys = sorted(StatsDict.keys())
	
					Percs = [0, 1, 2.5, 5, 10, 50, 90, 95, 97.5, 99, 100]
					# while we're at it, we get the scoreatpercentile for more than 
					# we will actually analyze in R.
					ScoreAtPercs = 	[scipy.stats.scoreatpercentile(StatsDict.get(x), perc) 		  for x in SortedKeys for perc in Percs]
					PercOfScore0s = [scipy.stats.percentileofscore(StatsDict.get(x), 0, kind='weak') for x in SortedKeys]
					# (results ordered saved in the order of SortedKeys)
	
					if Debugging:
						print len(Mins)
						print Mins
						print len(Maxs)
						print Maxs
						print len(Medians)
						print Medians
						print len(PercOfScore0s)
						print PercOfScore0s
						print len(ScoreAtPercs)			
						print ScoreAtPercs			
					# Stats complete. Print to file depending on how many records are left
	
					if (len(SpeciesRecords_groomed) != 0):   
						# we have data after grooming
						# do the full print
						PrintString = "%s\t%s\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (SpeciesPreviousLine, Len_all, Len_noDup, Len_groomed, Mins[0], Mins[1], Mins[2], Mins[3], Mins[4], Mins[5], Medians[0], Medians[1], Medians[2], Medians[3], Medians[4], Medians[5], Maxs[0], Maxs[1], Maxs[2], Maxs[3], Maxs[4], Maxs[5], PercOfScore0s[0], PercOfScore0s[1], PercOfScore0s[2], ScoreAtPercs[0], ScoreAtPercs[1], ScoreAtPercs[2], ScoreAtPercs[3], ScoreAtPercs[4], ScoreAtPercs[5], ScoreAtPercs[6], ScoreAtPercs[7], ScoreAtPercs[8], ScoreAtPercs[9], ScoreAtPercs[10], ScoreAtPercs[11], ScoreAtPercs[12], ScoreAtPercs[13], ScoreAtPercs[14], ScoreAtPercs[15], ScoreAtPercs[16], ScoreAtPercs[17], ScoreAtPercs[18], ScoreAtPercs[19], ScoreAtPercs[20], ScoreAtPercs[21], ScoreAtPercs[22], ScoreAtPercs[23], ScoreAtPercs[24], ScoreAtPercs[25], ScoreAtPercs[26], ScoreAtPercs[27], ScoreAtPercs[28], ScoreAtPercs[29], ScoreAtPercs[30], ScoreAtPercs[31], ScoreAtPercs[32])
	

				else:	
					# we don't have data after grooming, 
					# thus we calculate the stats only for _all and _noDup
					if (len(SpeciesRecords_noDup) != 0): 
						# just to be sure that we don't get
						# an error when a species only has NA or so.
						# even though noDup and all should always both exist
						# Extract latitude, longitude, and tmin from records
						Lats_all	 = [ x[1][0] for x in SpeciesRecords_all     ]
						Longs_all	 = [ x[1][1] for x in SpeciesRecords_all     ]
						Tmins_all	 = [ x[2] 	 for x in SpeciesRecords_all     ]
						Lats_noDup	 = [ x[1][0] for x in SpeciesRecords_noDup   ]
						Longs_noDup  = [ x[1][1] for x in SpeciesRecords_noDup   ]
						Tmins_noDup  = [ x[2] 	 for x in SpeciesRecords_noDup   ]
						
						# Assign extracted lat lon values to dictionary 
						StatsDict = {
							'Lats_all'	  	: Lats_all, 
							'Lats_noDup' 	: Lats_noDup,
							'Longs_all'	  	: Longs_all,
							'Longs_noDup'	: Longs_noDup
						}
						SortedKeys = sorted(StatsDict.keys())
		
						# Do the stats
						Mins    = [min(   StatsDict.get(x)) for x in SortedKeys]
						Medians = [median(StatsDict.get(x)) for x in SortedKeys]
						Maxs    = [max(   StatsDict.get(x)) for x in SortedKeys]
						# (results ordered saved in the order of SortedKeys)
						
						# Assign extracted tmin values to dictionary 
						StatsDict = {
							'Tmins_all'	  	: Tmins_all,
							'Tmins_noDup'	: Tmins_noDup
						}
						SortedKeys = sorted(StatsDict.keys())
		
						Percs = [0, 1, 2.5, 5, 10, 50, 90, 95, 97.5, 99, 100]
						ScoreAtPercs = 	[scipy.stats.scoreatpercentile(StatsDict.get(x), perc) 		  for x in SortedKeys for perc in Percs]
						PercOfScore0s = [scipy.stats.percentileofscore(StatsDict.get(x), 0, kind='weak') for x in SortedKeys]
						# (results ordered saved in the order of SortedKeys)
						PrintString = "%s\t%s\t%s\t%s\t%.8f\t\t%.8f\t%.8f\t\t%.8f\t%.8f\t\t%.8f\t%.8f\t\t%.8f\t%.8f\t\t%.8f\t%.8f\t\t%.8f\t%.3f\t\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\t\t\t\t\t\t\t\t\t\t\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (SpeciesPreviousLine, Len_all, Len_noDup, Len_groomed, Mins[0], Mins[1], Mins[2], Mins[3], Medians[0], Medians[1], Medians[2], Medians[3], Maxs[0], Maxs[1], Maxs[2], Maxs[3], PercOfScore0s[0], PercOfScore0s[1], ScoreAtPercs[0], ScoreAtPercs[1], ScoreAtPercs[2], ScoreAtPercs[3], ScoreAtPercs[4], ScoreAtPercs[5], ScoreAtPercs[6], ScoreAtPercs[7], ScoreAtPercs[8], ScoreAtPercs[9], ScoreAtPercs[10], ScoreAtPercs[11], ScoreAtPercs[12], ScoreAtPercs[13], ScoreAtPercs[14], ScoreAtPercs[15], ScoreAtPercs[16], ScoreAtPercs[17], ScoreAtPercs[18], ScoreAtPercs[19], ScoreAtPercs[20], ScoreAtPercs[21])

					else: 
						print "Oh No! the species below had no data!"
						print SpeciesPreviousLine
						PrintString = "NA\n"

				# write the result of parsing to the summary file
				SpecSummFile.write(PrintString)

 				# * RESET CURRENT Tuple and append current line
				SpeciesLineList = []
				SpeciesLineList.append(LineTuple)
		
		else: # Do this block if 'NA' present in LineElements 
			BrokenLinesFile.write(Line + '\n')
			BrokenLinesNumber += 1

		# end of processing line: remember species processed.
		SpeciesPreviousLine = SpeciesThisLine


	else: 	 
		# exception to handle header and first line of first species

		# if we are reading the header, write it to file		
		Line         = Line.strip('\n')
		LineElements = Line.split(',')
		if LineElements[0] == '"Scientific_name_interpreted"':   # it is the header
			# while we're at it, write header of summary file
			SpecSummFile.write(SpecSummFileHeader)
			
			
		else:   
			# it is NOT the header, thus the first species line; treat as normal line
			# read line and check for NA. 
			# if NA present, write error, then try again next line.
			if not('NA' in LineElements):   # we encountered a good line
				SpeciesNumber += 1
				SpeciesThisLine = LineElements[0]
				
				#  get coordinate precision. This is the same as above.
				LatPrecision = 0
				LongPrecision = 0
				MinImprecision = ['166666', '333333', '666666', '833333']
				LatLength  = str(LineElements[2])
				LatLength  = LatLength.split('.')
				if (len(LatLength) == 2):
					if LatLength[1][1:6] in MinImprecision:
						LatPrecision = int(1)
					else: 									 
						LatPrecision = int(len(LatLength[1]))
				else:
					LatPrecision = int(0)
				LongLength  = str(LineElements[1])
				LongLength  = LongLength.split('.')
				if (len(LongLength) == 2):
					if LongLength[1][1:6] in MinImprecision:
						LongPrecision = int(1)
					else: 
						LongPrecision = int(len(LongLength[1]))
				else:
					LongPrecision = int(0)
				
				# extract data from line and make it into sensible tuple
				tmp = [float(LineElements[x]) for x in range(1,4)]
				LineTuple = (SpeciesThisLine, (tmp[1], tmp[0]), tmp[2], LatPrecision, LongPrecision, )

				# start list of tuples, where each tuple is a line from the file
				SpeciesLineList.append(LineTuple)
			
				# switch to BeyondStart after the first good line
				BeyondStart = True

				# end of processing good line: remember species processed.
				SpeciesPreviousLine = SpeciesThisLine
		
			else: # first line contains NA: store it, try again, 
				BrokenLinesFile.write(Line + '\n')
				BrokenLinesNumber += 1
				BeyondStart = False
	
	# score the progress.			
	PrintNumber += 1
	if PrintNumber == 1000:
		PrintNumber = 0
		TimeNow = time.time()
		print "Lines read: %8d;  Species encountered: %7d;  Current species: %s;\tTime elapsed: %f" \
		  		% (LineNumber, SpeciesNumber, SpeciesThisLine, TimeNow-TimeStart)			

# done. Close files.

InFile.close()				# csv input
SpecSummFile.close()		# output: species summary stats
BrokenLinesFile.close()		# output: broken lines




