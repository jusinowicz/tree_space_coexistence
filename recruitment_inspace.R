###############################################################################
#
# This is code to take data on seedling records in plots with x,y coordinate
# data and turn it into a data frame of spatially explicit seedling 
# recruitment. These routines only generate the data to be used with 
# file "" which will be used to analyze spatiotemporal coexistence 
# mechanisms
#
############

#=============================================================================
# Load the data for BCI
#=============================================================================
seedlingXY = read.csv('SEEDLING_PLOT_XY.csv', header=TRUE)
seedling.full=read.csv('BCI_SEEDLING.csv', header=T)

#============================================================================
#First, remove records that are questionable or represent dead/damaged 
#individuals. Note: these codes are based on the metadata. 
#============================================================================

seedling=subset(seedling.full, seedling.full$NOTE!= 'DS' | seedling.full$NOTE!= 'DE' | seedling.full$NOTE != 'AQ' | seedling.full$NOTE != 'FL')

#============================================================================
# Make a raw "seedling" table with only first appearances of seedlings. 
# Use this with the below code to get annual recruitment. 
#============================================================================

#Seedling IDs
tags=unique(seedling$TAG)
ntags=length(tags)

seedling.rec = matrix ( data=0, 1, 9) #Make the basic matrix
colnames(seedling.rec)=colnames(seedling) #Label columns 
index=1;
#Outermost loop to go through all the unique tags:
for (i in 1:ntags){
		seed.temp=subset(seedling, seedling$TAG == tags[i]) #Find the record
		plots=unique(seed.temp$PLOT)
		nplots=length(plots)
	#Loop through all of the plots that a tag appeared in, and find its 
	#earliest appearance in each plot (i.e. recruitment date). 
	#This is necessary because tags were recycled after individuals died. 	
	for (j in 1:nplots){
			seed.plot.temp=subset(seed.temp, seed.temp$PLOT == plots[j])
			seedling.rec=rbind(seedling.rec, subset(seed.plot.temp, seed.plot.temp$YEAR == min(seed.plot.temp$YEAR)))
			}
}
#Clean up spurious records: 
seedling.rec=subset(seedling.rec, seedling.rec$YEAR != 0)

#To export this
write.table(seedling.rec, "seedling_recruits.csv", row.names=FALSE, col.names=TRUE)

#============================================================================
# This section takes the data object created in the previous section and 
# appends the X, Y data. 
#============================================================================

#There are 4 important factors here: Species, year, trap ID, and plot ID (each
#trap is associated with 3 seedling census plots)

species=unique(seedling.rec$SPECIES)
nspp=length(species)

years=unique(seedling.rec$YEAR)
nyrs=length(years)

traps=unique(seedling.rec$TRAP)
ntrap=length(traps)

plots=unique(seedlingXY$plot)
nplots=length(plots)

sppspace = matrix (data=0, nyrs*nplots, nspp+4)#create the main data object

index=1
#Outermost loop through years:
for(k in 1:nyrs){
	#Subset data by year.
	sppyear=subset(seedling.rec, seedling.rec$YEAR == years[k] ) 

	#Second loop through plots:
	for (j in 1:nplots) {
		#Subset by plot ID
		sppplot=subset(sppyear, sppyear$PLOT==plots[j] )

		#Third loop through species: 
		for (i in 1:nspp){
			#Subet by species
			spp=subset(sppplot,sppplot$SPECIES==species[i])
			#Now assign values based on where we are in each loop
			sppspace[index,1] = k
			sppspace[index,2] = j
			sppspace[index,3] = subset(seedlingXY, seedlingXY$plot == plots[j])$x
			sppspace[index,4] = subset(seedlingXY, seedlingXY$plot == plots[j])$y
			sppspace[index,4+i] = nrow(spp)
		}
	index=index+1 #Increment. This will represent total number of records
	}
}

#Final object:
colnames(sppspace)=c('YEAR', 'PLOT', 'X','Y',as.character(species))
write.table(sppspace, "seedling_rec_inspace.csv", row.names=FALSE, col.names=TRUE)
