#!/bin/bash/python

#==========================================
# QUANTIFICATION OF TURNOVER ALONG GRADIENTS
# By Melissa Chen
# Version 19-April-2020
# Built under python 3.7.14
# Created April 10th, 2017
#==========================================

##### TO DELETE####
# Not checked yet:
# changed taxaTable and metadata to taxaTablePWD and metadataPWD
# added minCountTable and minCountBin
# added minCountOTU filtering in loading of table
# added manual sigthresh
# got rid of tree
# fixed the way ubiq was calculated; fixed math but also added in a requirement that all ubiq must have half of observations as NON zeros.
# changed flow of inter blooms; was impossibly to reach before. Now, it is assessed with other bloom types
# changed order of saving typeoutput; was extraneous before.
# added bestfitTies file which tells you about warnings
# added levene's test to test difference in variances; also introduced new critpVar** need to assess side effects
# Maybe get rid of boundaries.txt in future because it seemed extraneous; will have to edit plotting script
# Got rid of checkValue because it didn't work as it should; not sure how it should work
# Got rid of suppress printing bins because it didn't seem to be necessary; might as well print everything.
# changed thresh test to <= for typeTaxa so that you can set to '0' if you want it to be TRULY absent
# tried using man whitney U test for means comparisons instead, since the data is very non-parametric. Seems to have more power this way; less things are classified as "noclass"
# sometimes manwhitney U test is too sensitive, esp during zero-inflated groups. So changed the requirements so that difference of means must be sig AND one of them must be greater than thresh
# Changed OTUTabletext.txt into the "after" OTU table, not before.Delete original text one
# Edited loading of OTU table so you *could* use non-qiime formats
# 11aug: added minSamplePres to get rid of OTUs that are exclusively in X number of samples (default at LEAST 3)
######

# This script takes an OTU table (output from QIIME) and mappingfile with salinity (QIIME-compatible) and produces the following text outputs:
	# 1. modelBoundaries_type.txt (Rows are unique OTUs; headers are 'OTU ID', 'taxonomy', 'type', 'typeSimple', 'X', 'Y', 'sigAB', 'sigBC', 'sigAC', 'boundaries', 'boundariestwo', 'bloom' )
	# 2. boundaries.txt (a list of all boundaries)
	# 3. taxa_abundances_across_salinity.txt (A table: Rows are taxa and headers are salinity; shows abundance across salinity (averaged for replicate salinities))
	# 4. types_across_salinity[_all/condensed].txt (A table: Rows are 'types' and headers are salinity; shows abundance of each type across salinity 
	# 5. OTUTablebyType.txt (OTU Table where rows are 'type' and columns are OTU IDs; meant count data)
	# 6. taxaIDLegend.txt (List of all unique OTUs and their assigned taxonomy)
	# 7. gradient.txt (List of names of gradient)
	# 8. LOG.txt (Printed output of settings for that run)
	
# What QTAG does is:
	# For each OTU in taxasummaries, finds a best-fit 3-piece model where each piece is the mean relative abundance for that region
		# - iterates through all X and Y (where X < Y, |X-Y| > diffXY, X > minVal and Y < maxVal)
	# Then, QTAG decides whether it is a Low,Inter, or High specialist (or none of them)
		# It does this by comparing the mean relative abundance of three groups (groupA, groupB, groupC),
		# The three groups are separated by boundaries X and Y.
		# Then, it uses Welch's T-test to see if there are significant differences in the means between groupA, groupB, and groupC
		# eg. if groupA > groupC (significant), then it is classified as Low.
		# It uses X and Y to calculate the places were each OTU seems to be 'turning over' (ie, changing in abundance significantly) 

# REQUIREMENTS:
	# Imports:
		# numpy
		# math
		# scipy stats
		# argparse
		# os
		# subprocess
		# time sleep
		#sys
	# Optional dependencies:
		# R script for plotting
		
		

#==========================================

import math 
import numpy # For standard deviation
	# Required for BINNING and BINNING SALINITY sections
from scipy import stats # For welch's t-test
	# Required TYPETAXA sections
# from graphics import *
import argparse # for user import
import os
import subprocess
from time import sleep
import sys
import copy # to copy lists and dictionaries
from progressbar import ProgressBar


#==========================================
# FUNCTIONS TO LOAD FILES

def makeTaxaSummaries(taxaTablePWD,metadataPWD):
	# Input is file path for OTU table and metadata
	# Output is 
		# (a) Dictionary where keys are OTU IDs and content is two lists. (taxaTable)
		# 		The first list contains the relative abundance of that OTU across salinity
		# 		The second list contains the corresponding salinities in the first list
		# (b) Dictionary where keys are OTU IDs and content is taxonomy (taxaIDs)
		# 		Taxonomies are taken from the observation metadata from the OTU table
	global metadata_name # This is the header name in the metadata file of the values to be assessed
	global minCountOTUinSample
	global minCountTable
	global minSamplePres
	print("Making and loading taxa table and metadata table...")
	taxaOpen = open(taxaTablePWD, 'r')
	taxaOpenTemp = []
	for i in taxaOpen:	# Read each line and concatenate into single file
		taxaOpenTemp += [i] # A single string; OTU table should be Unix(LF) format (\n at ends)
	taxaOpen.close()
	tempList =[] # to split lines
	for j in taxaOpenTemp: # Each line split by "\n"
		tempLine = j.strip()
		tempList += [tempLine.split('\t')] # Make a list of lists; each smaller list is abundance data for each OTU
	while '#OTU ID' not in tempList[0]:
		del tempList[0] # Deletes first row until the first row is #OTU ID
	# Sort information in relavent piles
	taxaIDs = {}
	taxaCountsTemp = {}
	first = True
	for y in tempList: # tempList is a list of lists; smaller list is OTUID,+ abundance data
		if y[0] == '#OTU ID': # The first line should be #OTU ID
			sites = y[1:len(y)-1] # So we will take the site names from the first row
		else: # Every other row is abundance data
			taxaIDs[y[0]] = y[len(y)-1] # Make dictionary of taxaIDs for later
			for x in range(len(y)): # Make file of 'total' counts for each OTU
				if (x != 0) and (x != (len(y)-1)): # If it's not the first element (which is the OTU ID) or last element (which is the taxonomy), then add to list
					if first: # First is made 'equal', but everything else is just appended.
						taxaCountsTemp[str(x)] = float(y[x])
					else:
						taxaCountsTemp[str(x)] += float(y[x])
			first = False
		# Output of this loop is taxaCountsTemp (a dictionary of abundance data), taxaIDs, and sites
	headers = sites[:]
	taxaTable = {}
	for i in range(len(tempList)):
		taxaTable[tempList[i][0]] = [[]] # Make empty list of lists for each OTU ID
		if tempList[i][0] == '#OTU ID': # 
			pass
		else:
			for j in range(1,len(tempList[i])-1):
				# Sum of all abundances to make relative abundances
	# 				sumAbund = int(taxaCountsTemp[str(j)])
				value = float(tempList[i][j])
				if value < float(minCountOTUinSample):
					value = 0.0	
				taxaTable[tempList[i][0]][0].append(value)
				# Save values as relative abundances instead of absolute ones
	# Now get rid of low abundance in table
	taxaTableFilt = deleteLowAbund(taxaTable,minCountTable,minSamplePres)
	# Get total counts
	totalCounts = [ 0 for i in range(len(sites))]
	for i in range(len(sites)):
		for taxa in taxaTableFilt:
			totalCounts[i] += taxaTableFilt[taxa][0][i]
	# Convert to relative abundance
	taxaTableFinal =copy.deepcopy(taxaTableFilt)
	for OTU in taxaTableFilt.keys():
		tempOTUlist = [0 for i in range(len(sites))]
		for i in range(len(sites)):
			tempOTUlist[i] = float(taxaTableFilt[OTU][0][i])/float(totalCounts[i])
		taxaTableFinal[OTU][0] = tempOTUlist		
	metadataOpen = open(metadataPWD, 'r')
	metadataOpenTemp = []
	for i in metadataOpen:
		metadataOpenTemp += [i]
	metadataOpen.close()
	tempMeta =[]
	for j in metadataOpenTemp:
		tempLine = j.strip()
		tempMeta += [tempLine.split('\t')]
	positionSal = tempMeta[0].index(metadata_name)
	metadata = []
	for line in tempMeta:
		metadata.append([line[0],line[positionSal]])
	# Now, change key names so they're not longer sites; they're gradient values
	for site in metadata:
		sites = [site[1] if x == site[0] else x for x in sites]
	# Make proper format, but with site names instead of numbers
	for taxa in taxaTableFinal:
		taxaTableFinal[taxa].append(sites)
	# Make abundance values integer as well
	for x in taxaTableFinal:
		taxaTableFinal[x] = [[float(strings) for strings in keys] for keys in taxaTableFinal[x]]
	printOTUTable(taxaTableFilt,headers,taxaIDs)
	return taxaTableFinal,taxaIDs

#==========================================
# Delete certain taxa based on total abundance; must be sufficiently abundant in order to work

def deleteLowAbund(taxaTable,minCountTable,minSamplePres): # Change to absolute threshold, or get rid of entirely. Get rid of show ups in < 3 samples
	newTaxaSummary = {}
	for taxa in taxaTable:
		nonzeroCount = 0
		presCount = 0
		for i in taxaTable[taxa][0]:
			nonzeroCount += i
			if i > 0:
				presCount += 1
		if int(nonzeroCount) >= int(minCountTable) and int(presCount) >= int(minSamplePres):
			newTaxaSummary[taxa] = taxaTable[taxa]
	return newTaxaSummary

#==========================================
# Msc functions

def average(listValues): # shortcut for 'average'
	# finds average for list of numbers
	if type(listValues) in [int,float]:
		listValues = [listValues]
	if len(listValues) == 0:
		return None
	else:
		average = float(sum(listValues))/len(listValues)
		return average
		
def is_numeric(s): # Is it a number
	try:
		float(s)
		return True
	except ValueError:
		return False
		
def printOTUTable(taxaTableFilt,headers,taxaIDs):
	OTUtabletoPrint = open('OTUTableText.txt','w')
	toPrint = '#OTUID'
	for head in headers:
		toPrint += '\t' + head
	toPrint += '\t'+'taxonomy'+'\n'
	for taxa in taxaTableFilt:
		toPrint += taxa
		for val in taxaTableFilt[taxa][0]:
			toPrint += '\t' + str(val)
		toPrint += '\t' + taxaIDs[taxa]
		toPrint += '\n'
	OTUtabletoPrint.write(toPrint)
	OTUtabletoPrint.close()
	print("Printed OTU table")
		
#==========================================
# MAKING BINS FOR EACH SET OF BOUNDARIES-- USED BELOW

def sortBins(X, Y, listAbund): # Sorts list of abundances into bins according to X and Y and their associated salinities
	lista = []
	listb = []
	listc = []
	for i in range(len(listAbund[1])):
		if listAbund[1][i] < X:  # if X = Y, value will be sorted into X
			lista.append(listAbund[0][i])
		if listAbund[1][i] >= X and listAbund[1][i] <= Y:
			listb.append(listAbund[0][i])
		if listAbund[1][i] > Y:
			listc.append(listAbund[0][i])
	binValues = {'lista':lista, 'listb':listb, 'listc':listc}
	return binValues # output is a dictionary with abundance observations sorted into lists according to what salinity they were found in

#==========================================

# MAKE AND TEST A MODEL USING ERROR SQUARED

def makeAModel(X,Y, binValues): # output is a dictionary with the abundances at each salinity, and then the 'model', which is just means
	aMean = average(binValues['lista'])
	bMean = average(binValues['listb'])
	cMean = average(binValues['listc'])
	# Make a combined 'piecewise' list that has each value and its corresponding mean
	abundance = binValues['lista'] + binValues['listb'] + binValues['listc']
	means = [aMean]*len(binValues['lista']) + [bMean]*len(binValues['listb']) + [cMean]*len(binValues['listc'])
	combined = {'abundance': abundance, 'means': means}
	return combined
	
def meanErrorSquared(combined): # Calculated the MES of data vs model
	# find the sum of squared difference for each value from the mean in that bin
	errorSquared = 0
	n = range(len(combined['abundance']))
	for i in n:
		error = combined['means'][i] - combined['abundance'][i]
		sqerror = error**2
		errorSquared += sqerror
	MES = float(errorSquared)/len(combined['abundance'])
	return MES

##  scalebyDiff no longer relavent since we're doing percentiles now
# def scaleByDiff(X,Y,meanA,meanB,meanC): # Scale X and Y into B given differences between means A,B,C
# 	if meanA > meanC:
# 		diffX = meanA-meanB
# 		diffY = meanB-meanC
# 		diffXY = Y-X
# 		if diffX <= 0:
# 			finalBoundary = Y
# 		elif diffY <= 0:
# 			finalBoundary = X
# 		else:
# 			scaleFactor = diffX/(diffX+diffY)
# 			diffscaled = diffXY*scaleFactor
# 			finalBoundary = Y-diffscaled
# 	if meanC > meanA:
# 		diffX = meanB-meanA
# 		diffY = meanC-meanB
# 		diffXY = Y-X
# 		if diffY <= 0:
# 			finalBoundary = X
# 		elif diffX <= 0:
# 			finalBoundary = Y
# 		else:
# 			scaleFactor = diffY/(diffX+diffY)
# 			diffscaled = diffXY*scaleFactor
# 			finalBoundary = X+diffscaled
# 	if meanA == meanC:
# 		print("WARNING:ERROR IN CALCULATING SCALEBYDIFF. Two means are exactly the same. Should not occur because it is classified as Hi or Lo")
# 	return finalBoundary
	
# def getBounds(listAbund,boundPerc):
# 	global Low,Inter,High
# 	# getting boundaries
# 	# Make listAbund only non-zeros
# 	tempAbund = [x for x in listAbund[0] if x>0]
# 	tempSal = [listAbund[1][i] for i,x in enumerate(listAbund[0]) if x>0]
# 	sumAbund = sum(tempAbund)
# 	percBound = sumAbund*boundPerc
# 	orderIndex = numpy.argsort(tempSal).tolist()
# 	orderedSal = [tempSal[i] for i in orderIndex]
# 	orderedAbund = [tempAbund[i] for i in orderIndex]
# 	salBound = None
# # 	currentIndex =None
# 	currentIndexuppr = None
# 	currentIndexlwr = None
# 	currentSum = 0
# 	halfBoundPerc=percBound/2
# 	currentIndexuppr = len(orderIndex)
# 	while currentSum < halfBoundPerc:
# 		currentIndexuppr -= 1
# 		currentSum += orderedAbund[currentIndexuppr]
# 	salBounduppr = orderedSal[currentIndexuppr-1]
# 	#reset for lwr bound
# 	currentSum = 0
# 	currentIndexlwr = -1
# 	while currentSum < halfBoundPerc:
# 		currentIndexlwr += 1
# 		currentSum += orderedAbund[currentIndexlwr]
# 	salBoundlwr = orderedSal[currentIndexlwr+1]
# # 	if salBoundlwr > salBounduppr:
# # 		temp =salBoundlwr
# # 		salBoundlwr = salBounduppr
# # 		salBounduppr = temp
# # 		print(listAbund)
# 	return([salBoundlwr, salBounduppr])
	
def getBounds(listAbund,boundPerc):
	# getting boundaries
	# Make listAbund only non-zeros
	tempAbund = [x for x in listAbund[0] if x>0]
	tempSal = [listAbund[1][i] for i,x in enumerate(listAbund[0]) if x>0]
	sumAbund = sum(tempAbund)
	percBound = sumAbund*boundPerc
	orderIndex = numpy.argsort(tempSal).tolist()
	orderedSal = [tempSal[i] for i in orderIndex]
	orderedAbund = [tempAbund[i] for i in orderIndex]
	# salBound = None
	currentIndexuppr = None
	currentIndexlwr = None
	# halfBoundPerc=percBound/2
	currentIndexuppr = len(orderIndex)-1
	currentIndexlwr = 0
	allAbund= [orderedAbund[i] for i in range(currentIndexlwr, currentIndexuppr+1)]
	currentSum = 0
	salBoundlwr = None
	salBounduppr = None
	for i in allAbund:
		currentSum += i
	notClose = currentSum > percBound
	while notClose and currentIndexlwr<currentIndexuppr: # make sure the bounds don't overlap
		changeInSum_lwr = orderedAbund[currentIndexlwr+1]
		changeInSum_upr = orderedAbund[currentIndexuppr-1]
		newSum_all = currentSum -changeInSum_lwr -changeInSum_upr
		newSum_lwr = currentSum - changeInSum_lwr
		newSum_uppr = currentSum - changeInSum_upr
		diff_all = abs(newSum_all-percBound)
		diff_lwr = abs(newSum_lwr-percBound)
		diff_uppr= abs(newSum_uppr-percBound)
		diff_current = abs(currentSum-percBound)
		if newSum_all > percBound: # if we are still  over the thresh, keep narrowing the boundaries.
			currentSum = newSum_all
			currentIndexlwr += 1
			currentIndexuppr -= 1
		elif diff_current < min([diff_lwr,diff_uppr]): # if moving BOTH boundaries puts us under thresh, check if the best model is actually the current one
			notClose = False
		elif newSum_lwr >= percBound: # we know changing both sides is too much; if we just change lower threshold, are we still safe?
			# If moving the lower threshold IS safe, let's check if moving just the upper is also safe
			if newSum_uppr >= percBound: # if this is true, then you can move either boundary; let's move the "smaller" jump.
				if changeInSum_lwr < changeInSum_upr: #is the smaller jump the lwr one?
					currentIndexlwr += 1
					currentSum = newSum_lwr
				else: # or is the smaller jumpt moving the higher one? if so, do below
					currentIndexuppr -= 1
					currentSum = newSum_uppr
			else: # If we reach here, then we know moving the lower is safe, but moving the upper is NOT safe
				currentIndexlwr += 1
				currentSum = newSum_lwr
		else: # Here, we know that moving the lower is NOT safe
			if newSum_uppr >= percBound: # changing lower is not safe, we check if moving the upper is safe; if yes, then move upper
				currentIndexuppr -=1
				currentSum = newSum_uppr
			else: # Here, we know that NEITHER are safe, so let's choose the lesser of two evils. It will put us past the threshold, but next round the best model will be the current one.
				if diff_lwr < diff_uppr: 
					currentIndexlwr += 1
					currentSum = newSum_lwr
				else: 
					currentIndexuppr -= 1
					currentSum = newSum_uppr
	salBoundlwr = orderedSal[currentIndexlwr]
	salBounduppr = orderedSal[currentIndexuppr]
	return([salBoundlwr, salBounduppr])

def typeTaxa(X, Y, listAbund): # Uses Welch's t-test and bins above to classify taxa as Hi, Lo, or Inter (or other)
	# This FIRST determines whether anything is intermediate:
	#	compares groups from best-fit model
	#	There are five intermediate sub-types: InterRestricted, InterPeak, InterPeakHiToler, InterPeakLoToler, Interbloom
	# Then, it tests whether it is Hi or Lo:
	# 	There are three intermediate sub-types each: Hi/LoRestricted, Hi/LoPeak, Hi/LoBloom
	# Finally, it tests whether the remaining groups are 'ubiquitous':
	#	Everything else is kept as 'noclass'
	global threshold
	global Low,Inter,High
	global ubiqOverlapThresh
	global propUbiqNonzero
	global critp
	global critpVar
	global boundPerc
	typeOutput = {'boundaries':[], 'type': '', 'typeSimple':'', 'meanA': [], 'meanB': [], 'meanC': [], 'X': [], 'Y': [], 'sigAB': [], 'sigBC': [], 'sigAC': [], 'bloom': 'No'} # Empty dictionary
	binValues = sortBins(X,Y,listAbund) # Use function above to create dictionary with abundance and salinity information
	typeOutput['boundaries'] = getBounds(listAbund,boundPerc) # get bounds out of the way
	# Find out whether someting is Hi-,Lo-, or Inter- specific
	groupA = binValues['lista']
	groupB = binValues['listb'] # should be at least 3 values in each list
	groupC = binValues['listc']
	meanA = average(groupA)
	meanB = average(groupB) 
	meanC = average(groupC)
	# Find threshold by using proportion of max, if necessariy
	if threshold[0] == True:
		maxAbund = max(groupA + groupB + groupC)
		thresh = maxAbund*threshold[1]
	else:
		thresh = threshold[1]
	# Calculate variance, but first test if each combination of groups has 0 variance.
	# If variance is 0, set p** to 1, which is maximum
	if average(groupA) == average(groupC) and numpy.var(groupA) == 0 and numpy.var(groupC) == 0:
		pAC = 1
	else:
		pAC = stats.mannwhitneyu(groupA,groupC)[1] # p-value of A vs C
	if average(groupA) == average(groupB) and numpy.var(groupA) == 0 and numpy.var(groupB) == 0:
		pAB = 1
	else:
		pAB = stats.mannwhitneyu(groupA,groupB)[1] # p-value of A vs B
	if average(groupB) == average(groupC) and numpy.var(groupB) == 0 and numpy.var(groupC) == 0:
		pBC = 1
	else:
		pBC = stats.mannwhitneyu(groupB,groupC)[1] # p-value of B vs C
	sigAB = pAB < critp # True if significant
	sigBC = pBC < critp # True if significant
	sigAC = pAC < critp # True if significant
	Stda = numpy.std(groupA) # For 'bloom' test
	Stdb = numpy.std(groupB)
	Stdc = numpy.std(groupC)
	typeOutput['sigAB'] = pAB
	typeOutput['sigBC'] = pBC
	typeOutput['sigAC'] = pAC
	typeOutput['meanA'] = meanA
	typeOutput['meanB'] = meanB
	typeOutput['meanC'] = meanC
	typeOutput['X'] = X
	typeOutput['Y'] = Y
	typeOutput['bloom'] =  'No'
	isInter = False	# See if there is an intermediate community or not; if there isn't, then I compare just X and C.
	isLow = False
	isHigh = False
	if sigAB and sigBC: # When the middle group is REAL (and not just 1 number), and it is significantly different than both flanking groups
		if meanB > meanA and meanB > meanC and meanB > thresh: # intermediate species
			isInter = True
# 			typeOutput['boundaries'] = [X,Y]
			typeOutput['typeSimple'] = Inter + 'Restricted'
			if meanA <= thresh and meanC <= thresh: # "very" intermediate; that is, the abundances are basically 0 on either side 
				typeOutput['type'] = Inter+'Restricted'
			elif meanA <= thresh and meanC > thresh: # leaning towards being hi
				typeOutput['type'] = Inter+'PeakHiToler'
			elif meanA > thresh and meanC <= thresh: # leaning towards being lo
				typeOutput['type'] = Inter+'PeakLoToler'
			else: # both meanA and meanC are larger than 0
				typeOutput['type'] = Inter+'Peak'
		elif meanB < meanA and meanB < meanC and meanA > thresh and meanC > thresh: # inv-inter water-- shouldn't exist according to hypothesis, but I put it in as a fail-safe
			isInter = True
			typeOutput['type'] = 'inv'+Inter
			typeOutput['typeSimple'] = 'noclass'
# 			typeOutput['boundaries'] = [X,Y]
		else:
 			pass	# All other combos mean the intermediates do NOT exist, so we can just compare X and C
# 		typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
	if isInter: # If the group is already classified as either Inter or inv-Inter, then we skip the next loop. If it is not classified, we continue.
		pass
	elif (meanA > meanC and sigAC and meanA > thresh): #or (meanB > meanC and sigBC and (meanB-meanC) > (meanB-meanA)): # More in fresh water and it's significant
		# Above, you can have EITHER X>C or Y>C but they must be significant, and if it's Y, the distance between Y and C must be greater than the distance between X and Y (to prevent Intermediate-looking ones)
		# Note that meanB-meanC should ALWAYS be greater than meanB-meanA because if meanB-meanA is negative, it means it's truly fresh!
		isLow = True
# 		typeOutput['boundaries'] = [scaleByDiff(X,Y,meanA,meanB,meanC)] # We scale by how significant each difference is. (eg. If a-b is very significant but b-c is not very significant, then the 'true' boundary is approximated to be closer to a-b than to b-c. 
# 		typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
		if meanC <= thresh:
			typeOutput['type'] = Low + 'Restricted'
			typeOutput['typeSimple'] = Low + 'Restricted'
		else:
			typeOutput['type'] = Low + 'Peak'
			typeOutput['typeSimple'] = Low + 'Restricted'
	elif (meanC > meanA and sigAC and meanC > thresh): #or (meanB > meanA and sigAB and (meanB-meanA) > (meanB-meanC)): # This is same process as above, except for marine samples
		isHigh = True
# 		typeOutput['boundaries'] = [scaleByDiff(X,Y,meanA,meanB,meanC)] # See above
# 		typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
		if meanA <= thresh:
			typeOutput['type'] = High + 'Restricted'
			typeOutput['typeSimple'] = High + 'Restricted'
		else:
			typeOutput['type'] = High + 'Peak'
			typeOutput['typeSimple'] = High + 'Restricted'
	if not any([isInter, isLow, isHigh]) : # if it has not been sorted yet:
		# Calculate whether the ranges of equal bins to overlap over the threshold set
		binThirds = sortBins(max(listAbund[1])*0.33,max(listAbund[1])*0.66,listAbund) 
		groupRanges = [[min(binThirds['lista']), max(binThirds['lista'])],[min(binThirds['listb']), max(binThirds['listb'])],[min(binThirds['listc']), max(binThirds['listc'])]]
		allOverlaps = []
		for i in range(0,len(groupRanges)):
			for j in range(0,len(groupRanges)):
				difflwr = max(0,groupRanges[j][0] - groupRanges[i][0])
				diffuppr = max(0,groupRanges[i][1] - groupRanges[j][1])
				totalDist = groupRanges[i][1]-groupRanges[i][0]
				if totalDist == 0:
					totalDist = 1
				finalOverlap = 1-(difflwr + diffuppr)/float(totalDist)
				allOverlaps.append(finalOverlap)
		# Also, calculate the number of zero and non-zero observations in ALL groups
		zerosprop = sum([1 for i in listAbund[0] if i == 0])/float(len(listAbund[0]))
		if zerosprop < propUbiqNonzero:
			ubiqZeroTestPass = True
		else:
			ubiqZeroTestPass = False
		if all([i >= ubiqOverlapThresh for i in allOverlaps]) and ubiqZeroTestPass: # passes both the overlap test and zero proportions test
			typeOutput['type'] = 'ubiquitous' # basically everywhere at similar levels
			typeOutput['typeSimple'] = 'noclass' # 'catch-all' for things that have NO significant differences between any of the groups
# 			typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
		else:
			##### TESTING OUT LEVENE'S TEST INSTEAD Of STD######
			#problem arises when two groups are EXACTLY the same; so let's check for that first.
			if len(set(groupA)) == 1 and len(set(groupB)) == 1 and set(groupA) == set(groupB):
				varABp = 1
			else:
				varABp = stats.levene(groupA,groupB)[1]
			if len(set(groupB)) == 1 and len(set(groupC)) == 1 and set(groupB) == set(groupC):
				varBCp = 1
			else:
				varBCp = stats.levene(groupB,groupC)[1]
			if len(set(groupA)) == 1 and len(set(groupC)) == 1 and set(groupA) == set(groupC):
				varACp = 1
			else:
				varACp = stats.levene(groupA,groupC)[1]
			sigvarAB = varABp < critpVar
			sigvarBC = varBCp < critpVar
			sigvarAC = varACp < critpVar
			if Stdb > Stda and Stdb > Stdc and sigvarAB and sigvarBC and meanC <= thresh and meanA <= thresh and meanB > meanC and meanB > meanA:
				isInter = True
				typeOutput['type'] = Inter+'Bloom'
				typeOutput['bloom'] = Inter
				typeOutput['typeSimple'] = Inter + 'Restricted'
# 				typeOutput['boundaries'] = [X,Y]
# 				typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
			elif Stdc > Stda and sigvarAC and meanC > meanA and meanC > meanB and meanA <= thresh:
				typeOutput['bloom'] = High
# 				typeOutput['boundaries'] = [scaleByDiff(X,Y,meanA,meanB,meanC)]
# 				typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
				typeOutput['type'] = High + 'Bloom'
				typeOutput['typeSimple'] = High + 'Restricted'
			elif Stda > Stdc and sigvarAC and meanA > meanC and meanA > meanB and meanC <= thresh:
				typeOutput['bloom'] =  Low
# 				typeOutput['boundaries'] = [scaleByDiff(X,Y,meanA,meanB,meanC)]
# 				typeOutput['boundaries'] = getBounds(listAbund,boundPerc)
				typeOutput['type'] = Low + 'Bloom'
				typeOutput['typeSimple'] = Low + 'Restricted'
			else:
				typeOutput['type'] = 'noclass'
				typeOutput['typeSimple'] = 'noclass' # 'catch-all' for things that have NO significant differences between any of the groups
# 				typeOutput['boundaries'] = ['',''] 
# 				typeOutput['boundaries'] = getBounds(listAbund,boundPerc) 
	return typeOutput # Output is a dictionary, in each taxa name there is 'type' and 'boundaries'; also, meanA, meanB, meanC, X, Y-- this is for downstream stuff

#==========================================
# GET TOLERANCE RANGES FOR ALL INDIVIDUALS

	
def getToleranceRange2(taxasummaries,stepDist,stepBinSize,taxaInfo):
	# 	intAbund_nozeros={'midSal': midSal, 'startSal': startSal, 'endSal': endSal} # dictionary of taxa; first is salinity midpoint; second is no zero averages; third is standardized by max; fourth is standardized by sum
	intAbund = {}
	for taxa in list(taxasummaries.keys()):
		# get boundaries to test
		minBound = math.floor(taxaInfo[taxa]['boundaries'][0])
		maxBound = math.ceil(taxaInfo[taxa]['boundaries'][1])
		# filter to ONLY within boundaries
		indexKeep = [i for i,x in enumerate(taxasummaries[taxa][1]) if x>=minBound and x<=maxBound]
		tempSal = [taxasummaries[taxa][1][i] for i in indexKeep]
		tempAbund = [taxasummaries[taxa][0][i] for i in indexKeep]
		# Store "extra" points outside of bounaries
		mscIndex = []
		for pair in enumerate(taxasummaries[taxa][1]):
			if pair[1]<minBound or pair[1]>maxBound:
				if taxasummaries[taxa][0][pair[0]] > 0:
					mscIndex.append(pair[0])
		mscPoints = [[x for i,x in enumerate(taxasummaries[taxa][1]) if i in mscIndex],[x for i,x in enumerate(taxasummaries[taxa][0]) if i in mscIndex]]
		# Need to create a sequence of bins
		startingSal = sequenceGenerator(minBound-stepBinSize+stepDist,maxBound,stepDist)
		# Get midpoint of segment we are averaging over
		midSal=[]
		for s in startingSal:
			tempmidSal = [minBound if i < minBound else i for i in sequenceGenerator(s, s+stepBinSize, stepDist)]
			midSal.append(average([maxBound if i>maxBound else i for i in tempmidSal]))
		# Get ranges for plotting later
		startSal = []
		endSal = []
		for i in range(len(midSal)):
			if i == 0:
				startSal.append(minBound)
			else:
				startSal.append(average([midSal[i-1],midSal[i]]))
			if i == (len(midSal)-1):
				endSal.append(maxBound)
			else:
				endSal.append(average([midSal[i],midSal[i+1]]))
		meanWithZero = []
		meanNoZero = []
		for s in startingSal:
				indexS = [i for i,x in enumerate(tempSal) if x >= s and x < (s + stepBinSize) ]
				tempAbund_current = [tempAbund[i] for i in indexS]
				meanWithZero.append(average(tempAbund_current))
				meanNoZero.append(average([ i for i in tempAbund_current if i>0]))
		# replace any Nones that arose from lack of samples
		withZero = [0 if v is None else v for v in meanWithZero]
		noZero = [0 if v is None else v for v in meanNoZero]
		peakWithZero = average([midSal[y] for y in [i for i, x in enumerate(withZero) if x == max(withZero)]])
		peakNoZero = average([midSal[y] for y in [i for i, x in enumerate(noZero) if x == max(noZero)]])
		# Save into dictionary
		intAbund[taxa] = {'midSal':midSal, 'startSal':startSal, 'endSal':endSal}
		intAbund[taxa].update({'abundances_withzero': withZero, 'abundances_withzero_maxstand':[x/max(withZero+mscPoints[1]) for x in withZero], 'abundances_withzero_sumstand':[x/sum(withZero+mscPoints[1]) for x in withZero]})
		intAbund[taxa].update({'abundances_nozero': noZero, 'abundances_nozero_maxstand':[x/max(noZero+mscPoints[1]) for x in noZero], 'abundances_nozero_sumstand':[x/sum(noZero+mscPoints[1]) for x in noZero]})
		# Get boundaries and midpoints for sorting later
		intAbund[taxa].update({'peakSal_withzero': peakWithZero, 'peakSal_nozero': peakNoZero})
		# add extreme values
		intAbund[taxa].update({'salPointOutliers':mscPoints[0], 'abundPointOutliers':mscPoints[1], 'abundPointOutliers_withzero_maxstand':[x/max(withZero+mscPoints[1]) for x in mscPoints[1]], 'abundPointOutliers_withzero_sumstand':[x/sum(withZero+mscPoints[1]) for x in mscPoints[1]]})
		intAbund[taxa].update({'abundPointOutliers_nozero_maxstand':[x/max(noZero+mscPoints[1]) for x in mscPoints[1]], 'abundPointOutliers_nozero_sumstand':[x/sum(noZero+mscPoints[1]) for x in mscPoints[1]]})
	return(intAbund)
	
#==========================================
# Take 'boundaries' that were significantly different and plot them to see how boundaries are distributed across salinities

def summaryBoundaryTypes(taxaInfo):
	global minX,maxY
	transitionsList = []
	transitions = []
	for taxa in taxaInfo:
		transitionsList = taxaInfo[taxa]['boundaries']
		if isinstance(transitionsList, list):
			for i in transitionsList:
				if i == '':
					pass
				else:
					transitions.append(i)
		else:
			if transitionsList == '':
				pass
			else:
				transitions.append(transitionsList)	
	return transitions # Output is summary of all boundaries in all taxa, excluding 'Nones'

#==========================================
# BINNING

def sequenceGenerator(min,max,binSize): # Make sequence from max value, min value, and binsize. List includes min and max bin values. Integers.
	if max < min:
		print("ERROR, max is less than min when trying to generate sequence in sequenceGenerator")
		return
	breadth = max-min
	nbins = math.ceil(float(breadth)/float(binSize))
	nbins = int(nbins) # Make into integer from float
	seq = []
	for i in range(nbins+1):
		current = float(min) + float(binSize)*i
		seq.append(current)
	return seq

#==========================================
# Determine composition (low, inter, high etc) within each gradient bin

class compositionAtSalinity():
	def __init__(self):
		self.lowGradRestr = 0
		self.lowGradPeak = 0
		self.interGradRestr = 0
		self.interGradPeakHiToler = 0
		self.interGradPeakLoToler = 0
		self.interGradPeak = 0
		self.invInterGrad = 0
		self.highGradRestr = 0
		self.highGradPeak = 0
		self.noclass = 0
		self.ubiquitous = 0
		self.lowGradBloom = 0
		self.highGradBloom = 0
		self.interGradBloom = 0
	def addValue(self,type,value,bloom): # Add 'raw' relative abundances from original taxasummaries
		if type == Low + 'Restricted':
			self.lowGradRestr += value
		if type == Low + 'Peak':
			self.lowGradPeak += value
		if type == Inter + 'Restricted':
			self.interGradRestr += value
		if type == Inter + 'PeakHiToler':
			self.interGradPeakHiToler += value
		if type == Inter + 'PeakLoToler':
			self.interGradPeakLoToler += value
		if type == Inter + 'Peak':
			self.interGradPeak += value
		if type == High + 'Restricted':
			self.highGradRestr += value
		if type == High + 'Peak':
			self.highGradPeak += value
		if type == 'noclass':
			self.noclass += value
		if type == Low + 'Bloom':
			self.lowGradBloom += value
		if type == Inter + 'Bloom':
			self.interGradBloom += value
		if type == High + 'Bloom':
			self.highGradBloom += value
		if type == 'ubiquitous':
			self.ubiquitous += value
		if type == 'inv'+Inter :
			self.invInterGrad += value
	def reportComp(self):
		total = self.lowGradRestr + self.lowGradPeak + self.interGradRestr + self.interGradPeakHiToler + self.interGradPeakLoToler + self.interGradPeak + self.highGradRestr + self.highGradPeak + self.noclass + self.ubiquitous + self.invInterGrad + self.lowGradBloom + self.highGradBloom + self.interGradBloom
		if total == 0: # Aka, there are no samples that exist in the data
			summary = 'Pass' 
		else: # By standardizing with "total", it should be a new 'average' relative abundance for each salinity bin.
			lowGradRestr = float(self.lowGradRestr)/total
			lowGradPeak = float(self.lowGradPeak)/total
			interGradRestr = float(self.interGradRestr)/total
			interGradPeakHiToler = float(self.interGradPeakHiToler)/total
			interGradPeakLoToler = float(self.interGradPeakLoToler)/total
			interGradPeak = float(self.interGradPeak)/total
			highGradRestr = float(self.highGradRestr)/total
			highGradPeak = float(self.highGradPeak)/total
			noclass = float(self.noclass)/total
			ubiquitous = float(self.ubiquitous)/total
			invInterGrad = float(self.invInterGrad)/total
			lowGradBloom = float(self.lowGradBloom)/total
			highGradBloom = float(self.highGradBloom)/total
			interGradBloom = float(self.interGradBloom)/total
			summary = [lowGradRestr,lowGradBloom,lowGradPeak,interGradPeakLoToler,interGradRestr,interGradBloom,interGradPeak,interGradPeakHiToler,highGradPeak,highGradRestr,highGradBloom,noclass,ubiquitous,invInterGrad]
		return summary # Output is list in order listed above.


def countListEvenLevels(listToCount,AllOptions):
	returnList = {}
	for i in AllOptions:
		returnList[i] = 0
		for j in listToCount:
			if j == i:
				returnList[i] += 1
			else:
				pass
	return returnList

#=========================================================



parser = argparse.ArgumentParser(
	description="Bins and classifies OTUs according to gradient specialization")
parser.add_argument(
	'-t',
	'--taxaTable',
	help = "Taxa table file, text format. Header for OTUs must be '#OTU ID', and number of lines before this don't matter.",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'File where first column are site names and another column is gradient values. Header required.',
	required = True)
parser.add_argument(
	'-M',
	'--metadata_name',
	help = 'Header name of column in metadata file that specifies the variable you want, if included',
	required = True)
parser.add_argument(
	'-o',
	'--output_dir',
	help = 'Output directory [default: GradientProcessing]',
	required = False,
	default = 'GradientProcessing')
parser.add_argument(
	'--gradient',
	help = 'Names for three points of gradient-- must have low, intermediate, and high name (Eg. Fresh, Brackish, Marine). Comma separated. [default: Low,Inter,High]',
	required = False,
	type = str,
	default = 'Low,Inter,High')
parser.add_argument(
	'--min_threshold_proportion',
	help = 'Proportional threshold for OTU abundances to be considered "non-zero". If both min_threshold_proportion and min_threshold_constant is supplied, constant will be used. [default 0.10]',
	required = False,
	type = float,
	default = 0.10)
parser.add_argument(
	'--min_threshold_constant',
	help = 'Constant threshold for OTU abundances to be considered "non-zero". If both min_threshold_proportion and min_threshold_constant is supplied, constant will be used. [default None]',
	required = False,
	default = None)
parser.add_argument(
	'--minX',
	help = "Minimum value for boundary X. IMPORTANT: default rounds the minimum gradient value down to nearest whole number. There may be situations where this value should be set manually. [default is minVal, the minimum gradient value]",
	required = False,
	default = 'Check')
parser.add_argument(
	'--maxY',
	help = 'Maximum value for boundary Y IMPORTANT: default rounds the maximum gradient value up to nearest whole number. There may be situations where this value should be set manually. [default is maxVal,the maximum gradient value]',
	required = False,
	default = 'Check')
parser.add_argument(
	'--XYdiff',
	help = 'Difference in UNITS between X and Y. Note: this is not raw difference between X and Y. This means how many increments you want to differentiate X and Y. For example, if your increments are 0.1PSU and you want X and Y to be at least 2PSU apart, then this value is actually 2/0.1=20 [default 2 units]',
	required = False,
	type = float,
	default = 2)
parser.add_argument(
	'--division_Size',
	help = 'Division size for bins for graphing bar graph [default 2]',
	required = False,
	default = 2)
parser.add_argument(
	'--unit_Size',
	help = 'Unit size of gradient to iterate through when fitting model. [default 1]',
	required = False,
	default = 1)
parser.add_argument(
	'--ubiq_overlap_threshold',
	help = 'When testing for ubiquitous taxa, the percent overlap required between bins to be considered ubiquitous. [Default: 0.10]',
	required = False,
	default = 0.10)
parser.add_argument(
	'--prop_ubiq_nonzero',
	help = 'When testing for ubiquitous taxa, the proportion of samples that must have non-zero abundances for that taxa. Will influence whether taxa is ubiquitous or noclass [Default: 0.30]',
	required = False,
	default = 0.30)
parser.add_argument(
	'--minCountTable',
	help = 'Minimum number of reads found in entire OTU table for an OTU to be kept. [Default: 50]',
	required = False,
	default = 50)
parser.add_argument(
	'--minCountBin',
	help = 'Minimum number of OTUs in each bin (A, B, C) in order for that bin combination to be valid. (Equal or greater than minCountBin) [Default: 3]',
	required = False,
	default = 3)
parser.add_argument(
	'--minCountOTUinSample',
	help = 'Minimum number of reads of any OTU in a sample in order for that OTU to be considered "non-zero". If the number of reads of that OTU in a sample is LESS than this value, then it will be changed to have an abundance of zero. [Default: 5]',
	required = False,
	default = 5)
parser.add_argument(
	'--minSamplePres',
	help = 'Threshold for number of samples an OTU must appear in. [default 3]',
	required = False,
	default = 3)
parser.add_argument(
	'--critp',
	help = 'Significance threshold to use when comparing differences in bins: controls how stringent you want differences between groups to be. Default is 0.05/2 (boneferroni corrected)[Default: 0.025]',
	required = False,
	default = 0.025)
parser.add_argument(
	'--critpVar',
	help = 'Significance threshold to use when comparing differences in variance within bins: controls how stringent you want bloom classifications to be. Lower p-value will mean more noclass and less bloom types. Default is 0.05 [Default: 0.05]',
	required = False,
	default = 0.05)	
parser.add_argument(
	'--stepDist',
	help = 'When creating average relative abundances for each OTU, the step size to take when advancing bin position [Default: 1]',
	required = False,
	default = 1)	
parser.add_argument(
	'--stepBinSize',
	help = 'When creating average relative abundances for each OTU, the bin size over which to average relative abundances [Default: 5]',
	required = False,
	default = 5)
parser.add_argument(
	'--boundPerc',
	help = 'Cumulative percent of summed relative abundances to use as cut-off for boundary determination [Default: 0.95]',
	required = False,
	default = 0.95)		
parser.add_argument(
	'-R',
	'--R_script',
	help = 'If given, will graph outputs through R script. Need pathway to Rscript file. If not provided, will not graph.',
	required = False,
	default = False)
	
args = parser.parse_args()

taxaTablePWD = args.taxaTable
metadataPWD = args.metadata
metadata_name = args.metadata_name
gradient = args.gradient
minX = args.minX
maxY = args.maxY
XYdiff = float(args.XYdiff)
min_threshold_proportion = float(args.min_threshold_proportion)
min_threshold_constant = args.min_threshold_constant
divisionSize = float(args.division_Size)
unitSize = float(args.unit_Size)
output_dir = args.output_dir
# allBins = args.allBins_SUPPRESS
# condensedBins = args.condensedBins_SUPPRESS
# checkClass = args.check_class
# checkValue = args.check_value
ubiqOverlapThresh = float(args.ubiq_overlap_threshold)
propUbiqNonzero = float(args.prop_ubiq_nonzero)
minCountTable = int(args.minCountTable)
minCountBin = int(args.minCountBin)
minCountOTUinSample = int(args.minCountOTUinSample)
minSamplePres = args.minSamplePres
critp = float(args.critp)
critpVar = float(args.critpVar)
stepDist = float(args.stepDist)
stepBinSize = float(args.stepBinSize)
boundPerc = float(args.boundPerc)
R_script = args.R_script


#==========================================
# Step one: load files and set variables

gradient = gradient.split(',')
Low = gradient[0]
Inter = gradient[1]
High = gradient[2]

# unitSize is how big increments are; XYdiff is minimum distance between X and Y
XYdiffunit = float(XYdiff)*float(unitSize)

# Check if they want an even threshold or a proportional threshold for 'zero' observations when deciding between restricted or tolerant types
if min_threshold_constant == None:
	threshold = [True,min_threshold_proportion]
else:
	threshold = [False,min_threshold_constant]
	
# Make taxa summary using function
taxasummaries,taxaIDs = makeTaxaSummaries(taxaTablePWD, metadataPWD)	

maxVal = max(taxasummaries[list(taxasummaries.keys())[1]][1])
minVal = min(taxasummaries[list(taxasummaries.keys())[1]][1])

if str(minX) == 'Check':
	minX = math.floor(minVal)
else:
	minX = minX
	
if str(maxY) == 'Check':
	maxY = math.ceil(maxVal)
else:
	maxY = maxY
	
#==========================================
# Step two: Choose best fit model and assign specialist type

# Make a list for storing warnings that come up during script
warningsFile = {"MES_zero":{}, "Bestfit_ties":{}}

# Dictionary that will have taxa as keys and a short list of paired [X,Y]
bestfitXY = {}
taxaInfo = {}

print("Iterating through all combinations...")
# Go through each taxa and iterate through all combinations of boundaries X and Y
print("PROGRESS")
totalTaxa = len(taxasummaries.keys())
pbar = ProgressBar()
for taxa in pbar(taxasummaries):
	currentTaxa = list(taxasummaries.keys()).index(taxa)
	currentbest = 0 # going to compare 1/n of error squared because I know it can't be lower than 0. Conversely, not sure what maximum of n will be.
	modelDiffList = []
	listAbund = taxasummaries[taxa]
	for X in sequenceGenerator(minX,(maxY-XYdiffunit),unitSize):
		for Y in sequenceGenerator((X+XYdiffunit),maxY,unitSize):
			binValues = sortBins(X,Y,listAbund)
			if len(binValues['lista']) <= minCountBin or len(binValues['listb']) <= minCountBin or len(binValues['listc']) <= minCountBin: # if the bins have nothing in them, then don't use that bin combination
				pass
			else:
				combined = makeAModel(X,Y,binValues)
				MES = meanErrorSquared(combined)
				if MES == 0: # If the meanErrorSquared is EXACTLY 0, there is probably something wrong, OR it is in such low abundance it's not worth looking at
					print("WARNING: MES IS ZERO")
					warningsFile['MES_zero'][taxa] = [X,Y]
					pass
				else:
					invMES = 1/MES
					modelDiffXsq = (average(binValues['lista']) - average(binValues['listb']))**2
					modelDiffYsq = (average(binValues['listb']) - average(binValues['listc']))**2
					modelDiff = (modelDiffXsq + modelDiffYsq)
					if invMES > currentbest:
						bestfitXY[taxa] = [[X,Y]]
						currentbest = invMES
						modelDiffList = [modelDiff]
					elif invMES == currentbest:
						bestfitXY[taxa].append([X,Y])
						modelDiffList.append(modelDiff)
	if len(bestfitXY[taxa]) == 1:
		aveFirst = bestfitXY[taxa][0][0]
		aveSecond = bestfitXY[taxa][0][1]
	else:
		firstBoundary = []
		secondBoundary = []
		for i in bestfitXY[taxa]:
			firstBoundary.append(i[0])
			secondBoundary.append(i[1])
		maxDifferenceFoundPosition = [i for i in range(len(modelDiffList)) if modelDiffList[i] == max(modelDiffList)]
		differencesX = numpy.diff(firstBoundary)
		differencesY = numpy.diff(secondBoundary)
		differencesXY = list(differencesX) + list(differencesY)
		if len(maxDifferenceFoundPosition) > 1 and (True in [True for i in differencesXY if i > 1]): # If consecutive numbers AND same diff
			print("WARNING: SAME MODELDIFF FOR IDENTICAL MODELS")
			print(taxa, firstBoundary, secondBoundary)
			warningsFile['Bestfit_ties'][taxa] = [firstBoundary,secondBoundary]
		firstBoundaries = []
		secondBoundaries = []
		for i in maxDifferenceFoundPosition:
			firstBoundaries.append(firstBoundary[i])
			secondBoundaries.append(secondBoundary[i])
		aveFirst = average(firstBoundaries)
		aveSecond = average(secondBoundaries)
	taxaInfo[taxa] = typeTaxa(aveFirst,aveSecond,listAbund) # type is going to be list of type and boundaries
	
print("DONE ITERATIONS")
transitions = summaryBoundaryTypes(taxaInfo) # list of all boundaries
binRanges = sequenceGenerator(minVal,maxVal,divisionSize)


#==========================================
# Step three: extract other msc data 

#==========================================
# Get tolerance ranges
# tolranges = getToleranceRange(taxasummaries,stepDist,stepBinSize,minX,maxY)
tolranges = getToleranceRange2(taxasummaries,stepDist,stepBinSize,taxaInfo)

#==========================================
# Make list of percent composition of salinities

# Make a file out of the taxasummaries file where you sum up the relative abundances at each salinity for each organism
# maxValinity = max(binRanges)

compositionAtGradAll = [0]*(len(binRanges)-1) # Empty matrix
for i in range(len(binRanges)-1): # For each bin
	compositionAtGradAll[i] = compositionAtSalinity()
	for x in taxasummaries: # For each taxa
		type = taxaInfo[x]['type']
		bloom = taxaInfo[x]['bloom']
		for y in range(len(taxasummaries[x][1])): # Add the relative abundance to the class, "compositionAtSalinity"
			if taxasummaries[x][1][y] >= binRanges[i] and taxasummaries[x][1][y] < binRanges[i+1]:
				compositionAtGradAll[i].addValue(type,taxasummaries[x][0][y],bloom)
composition = [0]*len(compositionAtGradAll)
for i in range(len(compositionAtGradAll)):
	composition[i] = compositionAtGradAll[i].reportComp() # Report the raw numbers for composition at each salinity

compositionAtGradCondensed = [0]*(len(binRanges)-1) # Empty matrix
for i in range(len(binRanges)-1): # For each bin
	compositionAtGradCondensed[i] = compositionAtSalinity()
	for x in taxasummaries: # For each taxa
		type = taxaInfo[x]['typeSimple']
		bloom = taxaInfo[x]['bloom']
		for y in range(len(taxasummaries[x][1])): # Add the relative abundance to the class, "compositionAtSalinity"
			if taxasummaries[x][1][y] >= binRanges[i] and taxasummaries[x][1][y] < binRanges[i+1]:
				compositionAtGradCondensed[i].addValue(type,taxasummaries[x][0][y],bloom)
compositionCondensed = [0]*len(compositionAtGradCondensed)
for i in range(len(compositionAtGradCondensed)):
	compositionCondensed[i] = compositionAtGradCondensed[i].reportComp() # Report the raw numbers for composition at each salinity

			
#==========================================
# Making OTU presence/absence table


allTypes = []
for i in taxaInfo:
	allTypes.append(taxaInfo[i]['type'])
allTypesUnique = list(set(allTypes))

OTUtableByType = {}
uniqueOTUList = []
for t in allTypesUnique: 
	OTUtableByType[t] = []
	for OTU in taxaInfo.keys():
		uniqueOTUList.append(OTU)
		if taxaInfo[OTU]['type'] == t:
				OTUtableByType[t].append(OTU)
		else: 
			pass
for t in allTypesUnique:
	tempdict = countListEvenLevels(OTUtableByType[t], uniqueOTUList)
	OTUtableByType[t] = tempdict
	
#==========================================


#==========================================
# Step four: Output the numerical data so we can look at it manually.

# combined = result from makeAModel
# transitions = list of all boundaries; can do histogram with output
# bestfitXY = dictionary, where indices are 'taxa' and then it lists their boundaries
# listAbund = [abundance,salinity]
print("Printing and saving...")

os.system('mkdir ' + output_dir)
os.chdir(output_dir)

# write 'all' version
typeAbundanceAll = open('types_across_gradient_all.txt', 'w')
toWrite = 'Gradient' + '\t'
for i in range(len(binRanges)-1):
	toWrite += str(binRanges[i+1]) + '\t'
toWrite = toWrite.strip()
toWrite += '\n'
typeAbundanceAll.write(toWrite)
categories = [Low+"Restricted",Low+"Bloom",Low + "Peak",Inter+"PeakLoToler",Inter+"Restricted",Inter+"Bloom",Inter+"Peak",Inter+"PeakHiToler",High+"Peak",High+"Restricted",High+"Bloom",'noclass','ubiquitous','inv'+Inter]
for i in range(len(categories)):
	toWrite = categories[i] + '\t'
	for j in composition:
		if j == 'Pass':
			toWrite += '\t'
		else:
			toWrite += str(j[i]) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\n'
	typeAbundanceAll.write(toWrite)
typeAbundanceAll.close()

# write 'condensed' version
typeAbundance = open('types_across_gradient_condensed.txt', 'w')
toWrite = 'Gradient' + '\t'
for i in range(len(binRanges)-1):
	toWrite += str(binRanges[i+1]) + '\t'
toWrite = toWrite.strip()
toWrite += '\n'
typeAbundance.write(toWrite)
categories = [Low+"Restricted",Low+"Bloom",Low + "Peak",Inter+"PeakLoToler",Inter+"Restricted",Inter+"Bloom",Inter+"Peak",Inter+"PeakHiToler",High+"Peak",High+"Restricted",High+"Bloom",'noclass','ubiquitous','inv'+Inter]
positionCategories = [0,4,9,11]
for i in positionCategories:
	toWrite = categories[i] + '\t'
	for j in compositionCondensed:
		if j == 'Pass':
			toWrite += '\t'
		else:
			toWrite += str(j[i]) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\n'
	typeAbundance.write(toWrite)
typeAbundance.close()
	

taxaFile = open('taxa_abundances_across_gradient.txt', 'w')
firstLine = 'Gradient\t'
for x in list(taxasummaries.values())[0][1]:
	firstLine += str(x) + '\t'
firstLine = firstLine.strip()
firstLine += '\n'
taxaFile.write(firstLine)
for taxa in taxasummaries:
	toWrite = taxa + '\t'
	for i in taxasummaries[taxa][0]:
		toWrite += str(i) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\n'
	taxaFile.write(toWrite)
taxaFile.close()

listBoundaries = open('boundaries.txt','w')
for x in transitions:
	listBoundaries.write(str(x) + '\n')
listBoundaries.close()

bestfitBoundaries = open('modelBoundaries_type.txt', 'w')
header = 'taxa\ttype\ttypeSimple\tmeanA\tmeanB\tmeanC\tA\tB\tsigAB\tsigBC\tsigAC\tbloom\tboundaries\tboundariestwo\n'
bestfitBoundaries.write(header)
for taxa in taxaInfo:
	toWrite = taxa + '\t'
	toWrite += taxaInfo[taxa]['type'] + '\t'
	toWrite += taxaInfo[taxa]['typeSimple'] + '\t'
	toWrite += str(taxaInfo[taxa]['meanA']) + '\t'
	toWrite += str(taxaInfo[taxa]['meanB']) + '\t'
	toWrite += str(taxaInfo[taxa]['meanC']) + '\t'
	toWrite += str(taxaInfo[taxa]['X']) + '\t'
	toWrite += str(taxaInfo[taxa]['Y']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigAB']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigBC']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigAC']) + '\t'
	toWrite += str(taxaInfo[taxa]['bloom']) + '\t'
	for i in taxaInfo[taxa]['boundaries']:
		toWrite += str(i) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\n'
	bestfitBoundaries.write(toWrite)
bestfitBoundaries.close()

taxaIDsprint = open('taxaIDLegend.txt','w')
toWrite = []
for taxa in taxaIDs:
	toWrite = str(taxa) + '\t' + str(taxaIDs[taxa]) + '\n'
	taxaIDsprint.write(toWrite)
taxaIDsprint.close()

gradientTemp = open("gradient.txt", 'w')
gradientTemp.write(Low+','+Inter+','+High+','+metadata_name+'\n')
gradientTemp.close()

# Print OTU presence/absence table 
OTUtabletoPrint = open('OTUTablebyType.txt','w')
firstDone = False
for t in OTUtableByType:
	lineToPrint = ''
	if not firstDone:
		lineToPrint += "#OTUTable \t"
		for i in OTUtableByType[t].keys():
			lineToPrint += str(i) + "\t"
		lineToPrint = lineToPrint.strip()
		lineToPrint += "\n"
		firstDone = True
	lineToPrint += t +"\t"
	for i in OTUtableByType[t]:
		lineToPrint += str(OTUtableByType[t][i]) + "\t"
	lineToPrint = lineToPrint.strip()
	lineToPrint += "\n"
	OTUtabletoPrint.write(lineToPrint)
OTUtabletoPrint.close()

uniqueListOfOTUs = open('OTUs.txt','w')
toPrint = ''
# OTUs = set(uniqueOTUList)
for taxa in taxaInfo:
	toPrint += str(taxa) + '\n'
toPrint = toPrint.strip()
uniqueListOfOTUs.write(toPrint)
uniqueListOfOTUs.close()

os.system('mv ../OTUTableText.txt ./')

# Print tolerance ranges
tolrangesPrint = open('toleranceRanges.txt', 'w')
header = 'OTU\tmidSal\tstartSal\tendSal\tpeakSal_withzero\tpeakSal_nozero\tabundances_withzero\tabundances_withzero_maxstand\tabundances_withzero_sumstand\tabundances_nozero\tabundances_nozero_maxstand\tabundances_nozero_sumstand\tposition\n'
tolrangesPrint.write(header)
for taxa in tolranges:
	for i in range(len(tolranges[taxa]['midSal'])):
		toWrite = taxa + '\t'
		toWrite += str(tolranges[taxa]['midSal'][i]) + '\t'
		toWrite += str(tolranges[taxa]['startSal'][i]) + '\t'
		toWrite += str(tolranges[taxa]['endSal'][i]) + '\t'
		toWrite += str(tolranges[taxa]['peakSal_withzero']) + '\t'
		toWrite += str(tolranges[taxa]['peakSal_nozero']) + '\t'
		toWrite += str(tolranges[taxa]['abundances_withzero'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundances_withzero_maxstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundances_withzero_sumstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundances_nozero'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundances_nozero_maxstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundances_nozero_sumstand'][i]) + '\t'
		toWrite += 'withinBounds' + '\n'
		tolrangesPrint.write(toWrite)
	for i in range(len(tolranges[taxa]['salPointOutliers'])):
		toWrite = taxa + '\t'
		toWrite += str(tolranges[taxa]['salPointOutliers'][i]) + '\t' # for midsal
		toWrite += str(tolranges[taxa]['salPointOutliers'][i]) + '\t' # for startsal
		toWrite += str(tolranges[taxa]['salPointOutliers'][i]) + '\t' # for endSal
		toWrite += str(tolranges[taxa]['peakSal_withzero']) + '\t'
		toWrite += str(tolranges[taxa]['peakSal_nozero']) + '\t'
		toWrite += str(tolranges[taxa]['abundPointOutliers'][i]) + '\t' # raw with zero
		toWrite += str(tolranges[taxa]['abundPointOutliers_withzero_maxstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundPointOutliers_withzero_sumstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundPointOutliers'][i]) + '\t' # raw no zero
		toWrite += str(tolranges[taxa]['abundPointOutliers_nozero_maxstand'][i]) + '\t'
		toWrite += str(tolranges[taxa]['abundPointOutliers_nozero_sumstand'][i]) + '\t'
		toWrite += 'outlier' + '\n'
		tolrangesPrint.write(toWrite)
tolrangesPrint.close()

# Print warnings file
toPrint='OTU\tMessage\n'
for otu in warningsFile['Bestfit_ties']:
	toPrint += otu + '\t' + str(warningsFile['Bestfit_ties'][otu]) + '\t WARNING: SAME MODELDIFF FOR IDENTICAL MODELS \n'
for otu in warningsFile['MES_zero']:
	toPrint += otu + '\t' + str(warningsFile['MES_zero'][otu]) + '\t WARNING: MES is zero \n'
WARNINGS = open("warnings.txt",'w')
WARNINGS.write(toPrint)
WARNINGS.close()

# PRINTING LOG
LOG = open("LOG.txt", 'w')
LOG.write(output_dir + '\n')
# if notes == True:
# 	LOG.write(raw_input("Type short description or comment, then press enter:\n") + '\n')
LOG.write("Settings: \n"+ "minX ="+ str(minX) + "\nmaxY =" + str(maxY) + "\nXYdiff =" + str(XYdiff) +  "\nDivSize =" + str(divisionSize) + "\nunitSize = " + str(unitSize))
thresholdType = 'Proportion'
if not threshold[0] :
	thresholdType = 'Constant'
LOG.write("\nthreshold type: " + thresholdType + "\nthreshold: " + str(threshold[1]))
LOG.write("\nFull PWD: \ntaxasummaries = " + taxaTablePWD + "\nmetadata = " + metadataPWD )
LOG.write("\nLow: "+ Low)
LOG.write("\nInter: "+ Inter)
LOG.write("\nHigh: "+ High)
LOG.write("\nGradient Header: "+ metadata_name)
LOG.write("\nUbiquitous taxa overlap threshold: " + str(ubiqOverlapThresh))
LOG.write("\nProportion overlap needed between groups to be ubiquitous: " + str(propUbiqNonzero))
LOG.write("\nFiltering; number of OTUs in sample to not be re-set to zero: " + str(minCountOTUinSample))
LOG.write("\nFiltering; number of reads of an OTU per sample for that OTU to be kept: " + str(minCountTable))
LOG.write("\nMinimum number of samples in bin: " + str(minCountBin))
LOG.write("\nCritical p-value to compare means of two bins: " + str(critp))
LOG.write("\nCritical p-value to compare variances when assigning bloom types: " + str(critpVar))

LOG.close()


#==========================================
# Step five: Draw the histogram and the other graphs
# Using R
print("Drawing graphs...")

if isinstance(R_script,str):
	subprocess.call(['Rscript',str(R_script)])

print('DONE')

