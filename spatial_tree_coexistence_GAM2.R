###############################################################################
#
# This is code to analyze SPATIAL coexistence from spatiotemporal tree
# recruitment data. This code includes only the spatial coexistence mechanisms, 
# as per Peter Chesson( Chesson 2000, Snyder and Chesson 2004).
#
# Coexistence is calculated on a pairwise basis between all species. 
# Coexistence is measured in terms of the invasion growth rates, with
# the fitness component removed by first removing the mean of recruitment
# from the data itself. 
#
# This version krigs missing spatial data using gstat
#
# This code is grouped into [4] parts: 
#	1. Load the data and pick out useful species (based on criteria of total 
#		number of samples, and total number of years with samples).
#	2. Use a GAM to krige recruit counts in the missing spatial locations
#	3. Calcuate the spatial components for pairwise coexistence, without any 
# 	 	simulation of population dyanmics. These components include:
#		The resident stationary distribution in the absence of the invader. 
#		The spatial storage effect. 	
#		Non-linear competitive variance. 
#		Fitness-density covariance. 
#		Total invasion growth rates.
#	4. Calculate (and graph) overall distributions of pairwise coexistence, 
#		measured by AijAji. These are equivalent to the inverse of the 
# 		the IGR when all fitness differences have been removed (thus the 
#		the Aij represent niche differences).
#		A. This part also calculates each of the spatial coexistence mechanisms
#			as a proportionate contribution to the total IGR!  
#	
###############################################################################
#=============================================================================
# Load these libraries
#=============================================================================
library(MASS)
library(mgcv)
library(fields)
library(raster)
require(parallel)
#For multiple cores
nc <- detectCores()   ## cluster size, set for example portability
if (detectCores()>1) { ## no point otherwise
  cl <- makeCluster(nc) 
  ## could also use makeForkCluster, but read warnings first!
} else cl <- NULL



#=============================================================================
# Load the data for BCI (see recruitment_inspace.R for details on data format)
#=============================================================================
recruits =read.table( "seedling_rec_inspace.csv", header=T) 

#=============================================================================
# Get some important variables from full data
#=============================================================================
ndata = 4 #The number of data columns
nspp.full = dim(recruits)[2]-ndata
nyears.full = length (unique (recruits$YEAR))

#=============================================================================
# Filter the data according to sampling thresholds. E.g. total number of 
# samples, total number of years sampled. 
#=============================================================================

min.samp=25 #minimum samples
re.use = recruits [,which(colSums(recruits)>=min.samp)]

yr.samp=3 #minimum number of years
nspp.tmp= dim(re.use)[2]-ndata
spp.use = matrix(0, nspp.tmp,1) #store logical, whether species meets criteria
for (n in 1:nspp.tmp){
	#For each species, subset based on when a seedling was sampled
	re.tmp = subset(re.use,re.use[,(n+4)]>0)
	#See if those samples span multiple years
	if ( length(unique(re.tmp$YEAR))>yr.samp){spp.use[n]=1}	
}
spp.use.names = colnames(re.use)[c(matrix(1,ndata,1),spp.use) == 1]
#Use only those species that passed the multi-year test
re.use=re.use[,spp.use.names]
nspp = dim(re.use)[2] - ndata
nyears = nyears.full



#=============================================================================
# FUNCTION DEFINITIONS
#=============================================================================
#R equivalent for matlab's meshgrid
meshgrid=function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 

#Assign values to a location within a matrix based on the xy coordinates of 
#these values. This assumes that another array object exists with the x and 
#y grids. 
xy_to_lattice = function () {


}


#TGet the fitness-density covariance for a single
#2D lattice of sites. It is based on Snyder 2008, Snyder and Chesson 2004. 
#The function returns a matrix with the (average) covariance in both the x
#and y directions in the first row, and their standard deviations in the 
#second row. As long as these averages are withing some statistical interval
#of each other, we can assume that spatial processes are isotropic and the
#covariance is independent of spatial direction.
#
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	cc.spp 	The competition from the resident community
#	l1.spp 	The spatial average of intrinsic reproduction
#	D.spp 	The denominator of the derivative in the Taylor expansion
#	gr1.spp The invasion growth rate WITHOUT the cov(f,d)
#	sr.spp 	The survival of the invader. 
#	D sets whether it is 1D or 2D:

get.fd.cov2D = function (Fr.spp, cc.spp, l1.spp, D.spp, gr1.spp, sr.spp, a_rr, D=2){ 
	rd1 = dim(Fr.spp)[1]
	rd2 = dim(Fr.spp)[2]
	cv.temp.r = matrix(0,rd1,1)
	cv.temp.c = matrix(0,rd1,1)
	cv.temp = matrix(0,2,2)

	#Covariance across rows
	for (r in 1:rd1) { 

		rmFr1=mean(Fr.spp[r,],na.rm=T) 
		if (rmFr1 > 0 & !is.na(rmFr1)) {
			rmcc2=mean(cc.spp[r,],na.rm=T) 
			ei1=(Fr.spp[r,]/rmFr1-1)
			uijmuj=(cc.spp[r,]/rmcc2-1)
			uijmuj[uijmuj<0]=0
			e1=(ei1-uijmuj)
			#Need to remove 0s for FFT
			e1[is.na(e1)] = 0
			s_zeta=spectrum(e1,spans=3,taper=0, plot=F, na.action=na.omit)
			ls1=length(s_zeta$spec)
			#Rebuild dispersal kernel to accomodate width
			xx_pad=matrix(seq(-ls1,ls1))
			kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
			kd.tmp=kd.tmp/(sum(kd.tmp))
			fdc2=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(ls1)
			#Then multiply by average igr
			agr1=(gr1.spp-sr.spp)
			cv.temp.r[r] = Re(fdc2[1])*agr1} else {

			cv.temp.r[r]=NA
			}
			
	}
		cv.temp[1,1] = mean(cv.temp.r,na.rm=T)
		cv.temp[2,1] = sqrt(var(cv.temp.r,na.rm=T))


	#Covariance across columns
	for (r in 1:rd2) { 

		rmFr1=mean(Fr.spp[,r],na.rm=T)
		if (rmFr1 > 0 & !is.na(rmFr1) ) {
			rmcc2=mean(cc.spp[,r],na.rm=T) 
			ei1=(Fr.spp[,r]/rmFr1-1)
			uijmuj=(cc.spp[,r]/rmcc2-1)
			uijmuj[uijmuj<0]=0
			e1=(ei1-uijmuj)
			#Need to remove 0s for FFT
			e1[is.na(e1)] = 0
			s_zeta=spectrum(e1,spans=3,taper=0, plot=F, na.action=na.omit)
			ls1=length(s_zeta$spec)
			#Rebuild dispersal kernel to accomodate width
			xx_pad=matrix(seq(-ls1,ls1))
			kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
			kd.tmp=kd.tmp/(sum(kd.tmp))
			fdc2=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(ls1)
			#Then multiply by average igr
			agr1=(gr1.spp-sr.spp)
			cv.temp.c[r] = Re(fdc2[1])*agr1} else {

			cv.temp.c[r]=NA
			}
			
	}

	cv.temp[1,2] = mean(cv.temp.r,na.rm=T)
	cv.temp[2,2] = sqrt(var(cv.temp.r,na.rm=T))


	return(cv.temp)


}

#=============================================================================
# END FUNCTION DEFINITIONS
#=============================================================================



#=============================================================================
# MAIN BODY: Dynamics and calculations of Low-density Growth Rates (LGR)
#=============================================================================
#=============================================================================
# Invasion growth rates for pairwise case, using the spatial lottery model
# as the underlying model. With this model specified, this can all be done 
# analytically. 
#=============================================================================

#=============================================================================
# USER ASSIGNED VALUES
#=============================================================================
#Prefix for naming files
f.name1=c("ldg_BCI_pairwise_")

###Survival
sr=matrix(0.9, nspp,1) 

###Competition coeffcients 
#alphas=matrix( c(1,0.8,0.8, 0.8,1,0.8,0.8,0.8,1),3,3)

###Competition distance
#b_rr=1/(100*np) #Essentially global
b_rr=matrix(1,nspp,1)

###Dispersal distance
#a_rr=1/(100*np) #essentially global
a_rr=matrix(1,nspp,1)


#=============================================================================
# Internal variable calls
#=============================================================================


#Spatial dimensions
# Note: use the "digits" argument of "round" to help determine spatial extent
x.max = round(max(re.use$X),-2)
y.max = round(max(re.use$Y),-2)
x.min = round(min(re.use$X),-2)
y.min = round(min(re.use$Y),-2)

#sig.xy=0.1 #The significant digit of the spatial measurements: for BCI, 0.1
sig.xy=1

xx = seq(x.min, x.max, sig.xy)
yy = seq(y.min, y.max, sig.xy)
np.xx=length(xx)
np.yy=length(yy)

#Make the full array of spatial coordinates.
stc.temp=meshgrid(xx,yy)
#Convert the coordinates into an array for easier access
stc=array(c(matrix(stc.temp$x,np.xx,np.yy,byrow=T),matrix(stc.temp$y,np.xx,np.yy,byrow=T)),dim=c(np.xx,np.yy,2)) 

#Make a second coordinate system for all of the spatial kernels, with positive and negative
#values and the origin defined as the center of the plot. NOTE: This system is TWICE the size
#of the original. This is done to 0-pad in 2D for the FFTs later! 
xx.c = seq(-np.xx, np.xx, sig.xy)
yy.c = seq(-np.yy, np.yy, sig.xy)
xx.c=xx.c[1: (np.xx*2)]
yy.c=yy.c[1: (np.yy*2)]
stc.c.temp=meshgrid(xx.c,yy.c)
stc.c=array(c(matrix(stc.c.temp$x,np.xx*2,np.yy*2,byrow=T),matrix(stc.c.temp$y,np.xx*2,np.yy*2,byrow=T)),dim=c(np.xx*2,np.yy*2,2)) 

#=============================================================================
# KRIGING: These data actually sample very little of the overall space. Use 
# these data and the positions of adult stems to co-krig the missing 
# data. 
#=============================================================================

re.use.k=NULL
sxy = data.frame(expand.grid(xx,yy))
colnames(sxy)=(c('X','Y'))
re.use.k = sxy

for (n in 1:nspp) {
	print(n)
	sd.rec=data.frame(cbind(re.use[,1:4],re.use[,n+4]))
	colnames(sd.rec) = c(colnames(re.use)[1:4],"count")

	#With GAM
	rec.gam = bam(count~te(X,Y,k=20) ,data=sd.rec, cluster=cl)
	#New data
	rec.krig = as.vector(predict.gam(rec.gam, newdata=sxy))
	#rec.krig = data.frame(cbind(sxy, rec.krig))
	#colnames(rec.krig)[3]= c("count")
	re.use.k=cbind(re.use.k,rec.krig)
	colnames(re.use.k)[n+2] = colnames(re.use)[n+4]
	#To plot
	#Rasterize the extent
	#rec.rast = as.matrix(rasterFromXYZ(rec.krig))


}

re.use.k=data.frame(re.use.k)
#save(file="seedling_rec_kriged.var", "re.use.k")
load(file="seedling_rec_kriged.var")

#===================================================================================

#Components of LGR for each species pair, with kriged data,
#with dispersal limitation or competition kernels

l1.k=matrix(0,nspp,nspp)
y.k=matrix(0,nspp,nspp)
cr_mean.k=matrix(0,nspp,nspp)
D.k = matrix(0,nspp,nspp)
var_mu_Us.k=matrix(0,nspp,nspp)
cov_e_mu_Us.k=matrix(0,nspp,nspp)
cov_lam_vc.k=matrix(0,nspp,nspp)
Elam1.k=matrix(0,nspp,nspp) 
Elam2.k=matrix(0,nspp,nspp)
gr1.n.k=matrix(0,nspp,nspp)
gr1.k=matrix(0,nspp,nspp)


#Look for spatial anisotropy in the spatial fd_covariance:
anis.yes = matrix(0,nspp,nspp)

#=============================================================================
#	Calculate the PAIRWISE LDG  with dispersal LIMITATION
#	This loops through all possible species pairs and calculates each 
# 	component of the LDG. 
#=============================================================================



for ( s in 1:nspp) { 
	#Declare the 3 main spatial variables used here: 
		# Frs the intrinsic reproduction of both species (before competition)
		# cc the competition that the invader experiences from the resident
		# nr2_eq_all the analytical expression for resident equilibrium density
	Frs = array(c(matrix(0,np.xx,np.yy),matrix(0,np.xx,np.yy)),dim=c(np.xx,np.yy,2)) 

	#Make a 2D lattice of values for the invader
	#Create a matrix with columns x,y,z, where x,y are plot positions and 
	#z is the count at that position. 
	Fr1 = cbind(round(re.use.k[,1:2]), re.use.k[,s+2])
	colnames(Fr1) = c(colnames(re.use.k)[1:2],colnames(re.use.k)[s+2])
	#Use the rasterized extent to project the x,y coords to a raster
	Frs[,,1] = matrix(rasterFromXYZ(Fr1),np.xx,np.yy)
	#Standardize by the mean
	Frs[,,1] = Frs[,,1]/mean(Frs[,,1],na.rm=T)

	#Loop through residents
	for (rs in 1:(nspp)){

		print( paste( "Invader: ", s, "Resident: ", rs))

		#Make a 2D lattice of values for the resident
		#Create a matrix with columns x,y,z, where x,y are plot positions and 
		#z is the count at that position. 
		Fr2 = cbind(re.use.k[,1:2], re.use.k[,rs+2])
		colnames(Fr2) = c(colnames(re.use.k)[1:2],colnames(re.use.k)[rs+2])
		#Use the rasterized extent to project the x,y coords to a raster
		Frs[,,2]= matrix(rasterFromXYZ(Fr2),np.xx,np.yy)
		Frs[,,2] = Frs[,,2]/mean(Frs[,,2],na.rm=T)

		####Dispersal kernels and their Fourier transforms: 

		#RIGHT NOW, THESE ARE IDENTICAL
		kd=matrix(0,np.xx*2,np.yy*2)
		kd = a_rr[rs]/2*exp(-a_rr[rs]*(abs(stc.c[,,1])+abs(stc.c[,,2])))
		kd=kd/(sum(kd))
		fkd=fft(kd)#/(np+1)


		####Competition kernels and their Fourier transforms:
		#RIGHT NOW, THESE ARE IDENTICAL
		kc=matrix(0,np.xx*2,np.yy*2)
		kc = b_rr[rs]/2*exp(-b_rr[rs]*(abs(stc.c[,,1])+abs(stc.c[,,2])))
		kc=kc/(sum(kc))
		fkc=fft(kc)#/(np+1)


		#Need to make NAs 0 for the FFT: 
		Frs[,,2][is.na(Frs[,,2])] = 0

		# This is the solved version of each resident's equilibrium density 
		# when it is alone. 
		# This version is based on using the model in Snyder 2006 but solving 
		# for the equilbrium relative density following Snyder 2004.  
		e1=Frs[,,2]/mean(Frs[,,2],na.rm=T)-1 #Transformed environmental response
		mat.new = matrix(0,np.xx*2,np.yy*2) # Embed this in a lot of 0s for fft
		mat.new [floor(np.xx/2):floor((np.xx+np.xx/2)-1),floor(np.yy/2):
			floor((np.yy+np.yy/2)-1) ]	= e1	
		eft=fft(mat.new)
		nr2_ave = mean((Frs[,,2]),na.rm=T)/(1-sr[rs]) 

		# The equilibrium equation for the model with E in numerator only. 
		# i.e. nr(t+1) = sr * nr (t) + conv(kr, Er/ conv(Uir, nr) )
		nr2_tmp=Re(fft((1-sr[rs])*fkd*eft/(fkd*(1-sr[rs])*(fkc-1)-
			sr[rs]+1),inverse=T))/((np.xx*2+0.2*np.xx*2)*(np.yy*2+0.2*np.yy*2))

		# This is the solved version using the approach in Snyder 2006.
		# It includes essentially the same terms, but is solved as the 
		# stationary distribution of an iterated process. 
		nr2_tmp.f = 0			
		for( n in seq(10,1,-1)) { 
			nr2_tmp.f = nr2_tmp.f+((1-sr[rs])*fkd*(1-fkc)+sr[rs])^(n-1) *(1-sr[rs])*fkd*eft
		}
		nr2_tmp = Re( fft(nr2_tmp.f,inverse=T))/((np.xx*2+0.2*np.xx*2)*(np.yy*2+0.2*np.yy*2))
		
		# Now take out the piece that we actually need:
		nr2_eq_all =matrix(0,np.xx*2,np.yy*2)
		nr2_eq_all[floor(np.xx/2):floor((np.xx+np.xx/2)-1),floor(np.yy/2):
			floor((np.yy+np.yy/2)-1) ] = 
			nr2_tmp[floor(np.xx/2):floor((np.xx+np.xx/2)-1),floor(np.yy/2):
			floor((np.yy+np.yy/2)-1) ]
		
		#This value is now in the format of the scaled variable u = ni/mean(ni) -1
		#Convert it back to density units for the competition: 
		nr2_eq_all = ( nr2_eq_all +1 )# *nr2_ave

		#Competition based on the resident equilibrium
		cc.k=Re(fft(fft(nr2_eq_all)*fkc,inverse=T)/((np.xx*2+0.2*np.xx*2)*(np.yy*2+0.2*np.yy*2)))
		cc.k=rbind(cc.k[ceiling(np.xx):(np.xx*2), ],cc.k[1:floor(np.xx), ])
		cc.k=cbind(cc.k[ ,ceiling(np.yy):(np.yy*2) ],cc.k[ ,1:floor(np.yy)])

		# Now rescale cc. k, and add NAs back in so the 0s (which are not real 0s) 
		# don't influence calculations
		Frs[,,2]= matrix(rasterFromXYZ(Fr2),np.xx,np.yy)
		Frs[,,2] = Frs[,,2]/mean(Frs[,,2],na.rm=T)

		cc.k = cc.k [floor(np.xx/2):floor((np.xx+np.xx/2)-1),floor(np.yy/2):
			floor((np.yy+np.yy/2)-1) ]
		cc.k[is.na(Frs[,,2])]=NA

		#Equilibrium average competitive fitness
		l1.k[s,rs]=mean(Frs[,,1],na.rm=T)
		y.k[s,rs]=mean(nr2_eq_all,na.rm=T)
		cr_mean.k[s,rs]=mean(cc.k,na.rm=T)
		D.k[s,rs]= cr_mean.k[s,rs]
		
		#Non-linear competitive variance and the spatial storage effect
		#These terms represent the perturbation terms from Snyder 2008 
		var_mu_Us.k[s,rs] = var( c((cc.k)/mean(cc.k,na.rm=T)-1),na.rm=T )
		cov_e_mu_Us.k[s,rs] = cov( c(Frs[,,1]/mean(Frs[,,1],na.rm=T)-1), c( (cc.k)/mean((cc.k),na.rm=T)-1),use="complete.obs" )
		
		#Components of the IGR
		Elam1.k[s,rs]=l1.k[s,rs]*(1/D.k[s,rs]+(cr_mean.k[s,rs]*var_mu_Us.k[s,rs])/D.k[s,rs]^3
			-(cr_mean.k[s,rs]*cov_e_mu_Us.k[s,rs])/D.k[s,rs]^2)+sr[s]-1
		Elam2.k[s,rs]=(l1.k[s,rs]*(1/D.k[s,rs]) +sr[s]-1)^2+2*(l1.k[s,rs]*(1/D.k[s,rs]) +sr[s]-1)*
			(var_mu_Us.k[s,rs]-cov_e_mu_Us.k[s,rs])/D.k[s,rs]^4
		gr1.n.k[s,rs] = exp(Elam1.k[s,rs]-0.5*Elam2.k[s,rs])
		
		#The fitness-density covariance 
		tryCatch( { cov.tmp = get.fd.cov2D(Frs[,,1], cc.k, l1.k[s,rs], D.k[s,rs], gr1.n.k[s,rs], sr[s], a_rr[s])
								cov_lam_vc.k[s,rs]= cov.tmp	[1,1]	
								#Look for anisotropy: that is, look for a difference in cov 
								#between the x and y directions
								sd.tmp = max(cov.tmp[2,1],cov.tmp[2,2]) #Which SD is larger
								if ( abs(cov.tmp[1,2] - cov.tmp[1,1]) >= sd.tmp ) {anis.yes[s,rs]=1}
					} , error=function(e){} )
		
		#The full IGR
		gr1.k[s,rs] = gr1.n.k[s,rs]+cov_lam_vc.k[s,rs]
		
		#save(file="spatial_coexistence_disp_kriged_BCI3.var", "l1.k","y.k", "cr_mean.k","D.k","var_mu_Us.k","cov_e_mu_Us.k",
		#"cov_lam_vc.k","Elam1.k", "Elam2.k", "gr1.n.k", "gr1.k")
		save(file="spatial_coexistence_disp_kriged2_BCI1.var", "l1.k","y.k", "cr_mean.k","D.k","var_mu_Us.k","cov_e_mu_Us.k",
		"cov_lam_vc.k","Elam1.k", "Elam2.k", "gr1.n.k", "gr1.k")


	}
}


#Some after-analysis 
#Load the results from the kriged data set.
load("spatial_coexistence_global_kriged_BCI.var") #Global values .gk
load("spatial_coexistence_disp_kriged_BCI_norm1B.var") #Limited dispersal .k
#Load the results from the raw data.
load("spatial_coexistence_global_BCI3.var") #Global values .g
load("spatial_coexistence_disp_norm_BCI2.var") #Limited dispersal .p


#These are the most analagous to the AijAji in temporal paper. 

aijaji.gk = 1/gr1.gk*t(1/gr1.gk)
aijaji.p = 1/gr1.p*t(1/gr1.p)
aijaji.k = 1/gr1.k*t(1/gr1.k)
aijaji.g = 1/gr1.g*t(1/gr1.g)

aijaji.gk = aijaji.gk[lower.tri(aijaji.gk)]
aijaji.p = aijaji.p[lower.tri(aijaji.p)]
aijaji.k = aijaji.k[lower.tri(aijaji.k)]
aijaji.g = aijaji.g[lower.tri(aijaji.g)]

#Calculate the percentage contribution of each spatial mechanism, to each 
#pairwise value of AijAji. 

#=============================================================================
# Figure comparing histograms of AijAji and their median values, for each of the
# four different sets of assumptions: 1) Global dispersal/competition, Kriged data;
# 2) Limited dispersal/competition, Kriged data; 3) Global, Raw data; and 
# 4) Limited, Raw data. 
#=============================================================================

pdf(file='all4_aijaji_BCI1.pdf', height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=14)
par(mfrow=c(2,2))
hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2] , xlim=c(0,1.2), xlab=' ', main="Global, Kriged", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],
		breaks=60,plot=F)$counts)/length(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2]),digits=2)))

hist(aijaji.k[aijaji.k > 0 & aijaji.k < 1.2], xlim=c(0,1.2), xlab=' ', main="Limited, Kriged", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijaji.k[aijaji.k > 0 & aijaji.k < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijaji.k[aijaji.k > 0 & aijaji.k < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijaji.k[aijaji.k > 0 & aijaji.k < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijaji.k[aijaji.k > 0 & aijaji.k < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijaji.k[aijaji.k > 0 & aijaji.k < 1.2],
		breaks=60,plot=F)$counts)/length(aijaji.k[aijaji.k > 0 & aijaji.k < 1.2]),digits=2)))

hist(aijaji.g[aijaji.g > 0 & aijaji.g < 1.2] , xlim=c(0,1.2), xlab=' ', main="Global, Raw", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijaji.g[aijaji.g > 0 & aijaji.g < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijaji.g[aijaji.g > 0 & aijaji.g < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijaji.g[aijaji.g > 0 & aijaji.g < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijaji.g[aijaji.g > 0 & aijaji.g < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijaji.g[aijaji.g > 0 & aijaji.g < 1.2],
		breaks=60,plot=F)$counts)/length(aijaji.g[aijaji.g > 0 & aijaji.g < 1.2]),digits=2)))

hist(aijaji.p[aijaji.p > 0 & aijaji.p < 1.2], xlim=c(0,1.2), xlab=' ', main="Limited, Raw", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijaji.p[aijaji.p > 0 & aijaji.p < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijaji.p[aijaji.p > 0 & aijaji.p < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijaji.p[aijaji.p > 0 & aijaji.p < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijaji.p[aijaji.p > 0 & aijaji.p < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijaji.p[aijaji.p > 0 & aijaji.p < 1.2],
		breaks=60,plot=F)$counts)/length(aijaji.p[aijaji.p > 0 & aijaji.p < 1.2]),digits=2)))

dev.off()


#=============================================================================
# Figure comparing histograms of the amount that each spatial mechanism contributes 
# to AijAji in each pairwise interaction. Graphs shown for two of the four 
# scenarios: 1) Global, Kriged data; and 2) Limited, Kriged data.
#=============================================================================
var_mu_Us.gkS = l1.gk*cr_mean.gk/(D.gk)^3*var_mu_Us.gk 
var_mu_Us.kS = l1.k*cr_mean.k/(D.k)^3*var_mu_Us.k 
cov_e_mu_Us.gkS = l1.gk*cr_mean.gk/(D.gk)^2*cov_e_mu_Us.gk 
cov_e_mu_Us.kS = l1.k*cr_mean.k/(D.k)^2*cov_e_mu_Us.k 


par(mfrow=c(3,2))
hist(var_mu_Us.gkS[var_mu_Us.gk > 0 & var_mu_Us.gk < 3] , xlim=c(0,1), xlab=' ', main="NCV: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(var_mu_Us.kS[var_mu_Us.k > 0 & var_mu_Us.k < 3] , xlim=c(0,1), xlab=' ', main="NCV: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.gkS[cov_e_mu_Us.gk > -0.5 & cov_e_mu_Us.gk < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.kS[cov_e_mu_Us.k > -0.5 & cov_e_mu_Us.k < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.gk[cov_lam_vc.gk > -2 & cov_lam_vc.gk < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.k[cov_lam_vc.k > -2 & cov_lam_vc.k < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Limited, Kriged", 
	cex.lab=1.3, breaks=60)

hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2] , xlim=c(0,1.2), xlab=' ', main="NCV: Global, Kriged", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2],
		breaks=60,plot=F)$counts)/length(aijaji.gk[aijaji.gk > 0 & aijaji.gk < 1.2]),digits=2)))


#=============================================================================
# Figure comparing histograms of AijAji and their median values, for each of the
# four different sets of assumptions: 1) Global dispersal/competition, Kriged data;
# 2) Limited dispersal/competition, Kriged data; 3) Global, Raw data; and 
# 4) Limited, Raw data. 
#=============================================================================
#These are the most analagous to the AijAji in temporal paper. 
#Set the scaling constant D to be equal (reflecting equivalent fitness, competition, and
# equilibrium resident abundance). 
Ds = 1

Elam1w.k = l1.k*(1/Ds+(var_mu_Us.k)/Ds^3
			-(cov_e_mu_Us.k)/Ds^2)+0.99-1
#Elam2w.k = (var_mu_Us.k-cov_e_mu_Us.k)/Ds^4
Elam2w.k =  (l1.k*(1/Ds) +0.99-1)^2+2*(l1.k*(1/Ds) +0.99-1)*
			(var_mu_Us.k-cov_e_mu_Us.k)/Ds^4
gr1w.k = exp(Elam1w.k-0.5*Elam2w.k)+cov_lam_vc.k

Elam1w.gk = l1.gk*(1/Ds+(var_mu_Us.gk)/Ds^3
			-(cov_e_mu_Us.gk)/Ds^2)+0.99-1
#Elam2w.gk = (var_mu_Us.gk-cov_e_mu_Us.gk)/Ds^4
Elam2w.gk =  (l1.gk*(1/Ds) +0.99-1)^2+2*(l1.gk*(1/Ds) +0.99-1)*
			(var_mu_Us.gk-cov_e_mu_Us.gk)/Ds^4
gr1w.gk = exp(Elam1w.gk-0.5*Elam2w.gk)+cov_lam_vc.gk

Elam1w.p = l1.p*(1/Ds+(var_mu_Us.p)/Ds^3
			-(cov_e_mu_Us.p)/Ds^2)+0.99-1
Elam2w.p =  (l1.p*(1/Ds) +0.99-1)^2+2*(l1.p*(1/Ds) +0.99-1)*
			(var_mu_Us.p-cov_e_mu_Us.p)/Ds^4
gr1w.p = exp(Elam1w.p-0.5*Elam2w.p)+cov_lam_vc.p

Elam1w.g = l1.g*(1/Ds+(var_mu_Us.g)/Ds^3
			-(cov_e_mu_Us.g)/Ds^2)+0.99-1
Elam2w.g =  (l1.g*(1/Ds) +0.99-1)^2+2*(l1.g*(1/Ds) +0.99-1)*
			(var_mu_Us.g-cov_e_mu_Us.g)/Ds^4
gr1w.g = exp(Elam1w.g-0.5*Elam2w.g)+cov_lam_vc.g

aijajiw.gk = 1/gr1w.gk*t(1/gr1w.gk)
aijajiw.p = 1/gr1w.p*t(1/gr1w.p)
aijajiw.k = 1/gr1w.k*t(1/gr1w.k)
aijajiw.g = 1/gr1w.g*t(1/gr1w.g)

aijajiw.gk = aijajiw.gk[lower.tri(aijajiw.gk)]
aijajiw.p = aijajiw.p[lower.tri(aijajiw.p)]
aijajiw.k = aijajiw.k[lower.tri(aijajiw.k)]
aijajiw.g = aijajiw.g[lower.tri(aijajiw.g)]

#The figure: 
pdf(file='all4_aijaji_stand_BCI1.pdf', height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=14)
par(mfrow=c(2,2))
hist(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 1.2] , xlim=c(0,1.2), xlab=' ', main="Global, Kriged", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 1.2],
		breaks=60,plot=F)$counts)/length(aijajiw.gk[aijajiw.gk > 0 & aijajiw.gk < 1.2]),digits=2)))

hist(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 1.2], xlim=c(0,1.2), xlab=' ', main="Limited, Kriged", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 1.2],
		breaks=60,plot=F)$counts)/length(aijajiw.k[aijajiw.k > 0 & aijajiw.k < 1.2]),digits=2)))

hist(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 1.2] , xlim=c(0,1.2), xlab=' ', main="Global, Raw", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 1.2],
		breaks=60,plot=F)$counts)/length(aijajiw.g[aijajiw.g > 0 & aijajiw.g < 1.2]),digits=2)))

hist(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 1.2], xlim=c(0,1.2), xlab=' ', main="Limited, Raw", 
	xaxt='n', yaxt='n',cex.lab=1.3, breaks=60)
#Only put the median on the x axis
axis(side=1, labels=F, at=c(0,round(median(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 2]),2),1))
axis(side=1, labels=T, at=c(round(median(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 2]),2)))
axis(side=1, labels=T, at=c(0,1.2))
axis(side=2, labels=F, at=c(0,max(hist(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 1.2],breaks=60,plot=F)$counts)))
axis(side=2, at=max(hist(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 1.2],breaks=60,plot=F)$counts),
	labels=c(round(max(hist(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 1.2],
		breaks=60,plot=F)$counts)/length(aijajiw.p[aijajiw.p > 0 & aijajiw.p < 1.2]),digits=2)))

dev.off()
#=============================================================================
# Figure comparing histograms of the amount that each spatial mechanism contributes 
# to AijAji in each pairwise interaction expressed as a percentage. 
# Graphs shown for two of the four 
# scenarios: 1) Global, Kriged data; and 2) Limited, Kriged data.
#=============================================================================
#As a percentage of each pairwise aijaji:
#per_var.gk = 1 - (exp((l1.gk*(1/Ds+(-cov_e_mu_Us.gk)/Ds^2)+0.99-1) 
#			-0.5*(l1.gk*(1/Ds) +0.99-1)^2+2*(l1.gk*(1/Ds) +0.99-1)*
#			(-cov_e_mu_Us.gk)/Ds^4 ) + cov_lam_vc.gk )/gr1w.gk

per_var.gk = 1 - (((l1.gk*(1/Ds+(-cov_e_mu_Us.gk)/Ds^2)+0.99) 
			-0.5*(l1.gk*(1/Ds) +0.99-1)^2+2*(l1.gk*(1/Ds) +0.99)*
			(-cov_e_mu_Us.gk)/Ds^4 ) + cov_lam_vc.gk )/gr1w.gk

per_se.gk = 1 - ( ( (l1.gk*(1/Ds+(var_mu_Us.gk)/Ds^3)+0.99) 
			- 0.5*(l1.gk*(1/Ds) +0.99-1)^2+2*(l1.gk*(1/Ds) +0.99)*
			(var_mu_Us.gk)/Ds^4)+cov_lam_vc.gk )/gr1w.gk

per_fdcov.gk = 1 - (exp(Elam1w.gk-0.5*Elam2w.gk) ) / gr1w.gk


per_var.k = 1 - (((l1.k*(1/Ds+(-cov_e_mu_Us.k)/Ds^2)+0.99) 
			-0.5*(l1.k*(1/Ds) +0.99-1)^2+2*(l1.k*(1/Ds) +0.99)*
			(-cov_e_mu_Us.k)/Ds^4 ) + cov_lam_vc.k )/gr1w.k

per_se.k = 1 - ( ( (l1.k*(1/Ds+(var_mu_Us.k)/Ds^3)+0.99) 
			- 0.5*(l1.k*(1/Ds) +0.99-1)^2+2*(l1.k*(1/Ds) +0.99)*
			(var_mu_Us.k)/Ds^4)+cov_lam_vc.k )/gr1w.k

per_fdcov.k = 1 - (exp(Elam1w.k-0.5*Elam2w.k) ) / gr1w.k


per_var.gk = per_var.gk[lower.tri(per_var.gk)]
per_var.k = per_var.k[lower.tri(per_var.k)]
per_se.gk = per_se.gk[lower.tri(per_se.gk)]
per_se.k = per_se.k[lower.tri(per_se.k)]
per_fdcov.gk = per_fdcov.gk[lower.tri(per_fdcov.gk)]
per_fdcov.k = per_fdcov.k[lower.tri(per_fdcov.k)]



pdf(file='all4_coexistenceMech_stand_BCI1.pdf', height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=14)

par(mfrow=c(3,2))

hist(per_var.gk[per_var.gk >= 0 & per_var.gk <= 1] , xlim=c(0,1), xlab=' ', main="NCV: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(per_var.k[per_var.k >= 0 & per_var.k <= 1] , xlim=c(0,1), xlab=' ', main="NCV: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(per_se.gk[per_se.gk >= -1 & per_se.gk<=1] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(per_se.k[per_se.k >= -1 & per_se.k<=1] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(per_fdcov.gk[per_fdcov.gk >= 0 & per_fdcov.gk <=1] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(per_fdcov.k[per_fdcov.k >= 0 & per_fdcov.k <= 1] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Limited, Kriged", 
	cex.lab=1.3, breaks=60)



dev.off()

#=============================================================================
# Figure comparing histograms of the amount that each spatial mechanism contributes 
# to AijAji in each pairwise interaction. Graphs shown for two of the four 
# scenarios: 1) Global, Kriged data; and 2) Limited, Kriged data.
#=============================================================================
var_mu_Us.gkS = 1/(2)^3*var_mu_Us.gk 
var_mu_Us.kS = 1/(2)^3*var_mu_Us.k 
cov_e_mu_Us.gkS = 1/(2)^2*cov_e_mu_Us.gk 
cov_e_mu_Us.kS = 1/(2)^2*cov_e_mu_Us.k 


pdf(file='all4_coexistenceMech_stand_BCI1.pdf', height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=14)

par(mfrow=c(3,2))


pdf(file='all4_coexistenceMech_stand_BCI1.pdf', height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=14)

par(mfrow=c(3,2))
hist(var_mu_Us.gkS[var_mu_Us.gk > 0 & var_mu_Us.gk < 3] , xlim=c(0,1), xlab=' ', main="NCV: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(var_mu_Us.kS[var_mu_Us.k > 0 & var_mu_Us.k < 3] , xlim=c(0,1), xlab=' ', main="NCV: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.gkS[cov_e_mu_Us.gk > -0.5 & cov_e_mu_Us.gk < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.kS[cov_e_mu_Us.k > -0.5 & cov_e_mu_Us.k < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Limited, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.gk[cov_lam_vc.gk > -2 & cov_lam_vc.gk < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Global, Kriged", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.k[cov_lam_vc.k > -2 & cov_lam_vc.k < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Limited, Kriged", 
	cex.lab=1.3, breaks=60)

dev.off()


#=============================================================================
# Figure comparing histograms of the amount that each spatial mechanism contributes 
# to AijAji in each pairwise interaction. Graphs shown for two of the four 
# scenarios: 1) Global, Raw data; and 2) Limited,Raw data.
#=============================================================================
var_mu_Us.gS = 1/(2)^3*var_mu_Us.g
var_mu_Us.pS = 1/(2)^3*var_mu_Us.p
cov_e_mu_Us.gS = 1/(2)^2*cov_e_mu_Us.g 
cov_e_mu_Us.pS = 1/(2)^2*cov_e_mu_Us.p


par(mfrow=c(3,2))
hist(var_mu_Us.gS[var_mu_Us.g > 0 & var_mu_Us.g  < 3] , xlim=c(0,1), xlab=' ', main="NCV: Global, Raw", 
	cex.lab=1.3, breaks=60)
hist(var_mu_Us.pS[var_mu_Us.p > 0 & var_mu_Us.p < 3] , xlim=c(0,1), xlab=' ', main="NCV: Limited, Raw", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.gS[cov_e_mu_Us.g > -0.5 & cov_e_mu_Us.g < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Global, Raw", 
	cex.lab=1.3, breaks=60)
hist(cov_e_mu_Us.pS[cov_e_mu_Us.p > -0.5 & cov_e_mu_Us.p < 0.5] , xlim=c(-0.5,0.5), xlab=' ', main="SE: Limited, Raw", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.g[cov_lam_vc.g > -2 & cov_lam_vc.g < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Global, Raw", 
	cex.lab=1.3, breaks=60)
hist(cov_lam_vc.p[cov_lam_vc.p > -2 & cov_lam_vc.p < 2] , xlim=c(-0.5,0.5), xlab=' ', main="FDC: Limited, Raw", 
	cex.lab=1.3, breaks=60)

