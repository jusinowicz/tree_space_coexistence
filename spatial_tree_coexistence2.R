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
# 
#
# This code is grouped into [] parts: 
#	1. Load the data and pick out useful species (based on criteria of total 
#		number of samples, and total number of years with samples).
#	2. Calculate the resident distribution in the absence of the invader.
#		A. Create the spatiotemporal recruitment data for an arbitrary of 
#		generations based on repeated bootstrapping of the original data set.
#		B. Run the population dynamics for the resident. 
#
###############################################################################
#=============================================================================
# Load these libraries
#=============================================================================
library(MASS)
library(fields)
library(geoR)
library(raster)


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

g
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


#Components of LGR for each species pair, 
#without dispersal limitation or competition kernels

l1.p=matrix(0,nspp,nspp)
y.p=matrix(0,nspp,nspp)
cr_mean.p=matrix(0,nspp,nspp)
D.p = matrix(0,nspp,nspp)
var_mu_Us.p=matrix(0,nspp,nspp)
cov_e_mu_Us.p=matrix(0,nspp,nspp)
cov_lam_vc.p=matrix(0,nspp,nspp)
Elam1.p=matrix(0,nspp,nspp) 
Elam2.p=matrix(0,nspp,nspp)
gr1.n.p=matrix(0,nspp,nspp)
gr1.p=matrix(0,nspp,nspp)

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
	Fr1 = cbind(round(re.use[,3:4]), re.use[,s+ndata])
	colnames(Fr1) = c(colnames(re.use)[3:4],colnames(re.use)[s+ndata])
	#Define the extent of the spatial grid 
	fr.ext = extent(matrix(c(x.min,x.max,y.min,y.max),2,2,byrow=T))
	#Rasterize the extent
	fr.r = raster(fr.ext, ncol=np.yy, nrow=np.xx)
	#Use the rasterized extent to project the x,y coords to a raster
	Fr1.rast = rasterize(Fr1[,1:2],fr.r,Fr1[,3],fun=mean)
	#Extract the matrix of data points that correspond to sample locations
	Frs[,,1]= as.matrix(Fr1.rast)
	#Standardize by the mean
	Frs[,,1] = Frs[,,1]/mean(Frs[,,1],na.rm=T)

	

	#Loop through residents
	for (rs in 1:(nspp)){

		print( paste( "Invader: ", s, "Resident: ", rs))

		#Make a 2D lattice of values for the resident
		#Create a matrix with columns x,y,z, where x,y are plot positions and 
		#z is the count at that position. 
		Fr2 = cbind(re.use[,3:4], re.use[,rs+ndata])
		colnames(Fr2) = c(colnames(re.use)[3:4],colnames(re.use)[rs+ndata])
		#Use the rasterized extent to project the x,y coords to a raster
		Fr2.rast = rasterize(Fr2[,1:2],fr.r,Fr2[,3],fun=mean)
		#Extract the matrix of data points that correspond to sample locations
		Frs[,,2]= as.matrix(Fr2.rast)
		#Standardize by the mean
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
		#nr2_tmp=Re(fft((1-sr[rs])*fkd*eft/(fkd*(1-sr[rs])*(fkc-1)-
		#	sr[rs]+1),inverse=T))/((np.xx*2+0.2*np.xx*2)*(np.yy*2+0.2*np.yy*2))

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
		cc.p=Re(fft(fft(nr2_eq_all)*fkc,inverse=T)/((np.xx*2+0.2*np.xx*2)*(np.yy*2+0.2*np.yy*2)))
		cc.p=rbind(cc.p[ceiling(np.xx):(np.xx*2), ],cc.p[1:floor(np.xx), ])
		cc.p=cbind(cc.p[ ,ceiling(np.yy):(np.yy*2) ],cc.p[ ,1:floor(np.yy)])


		#Now add NAs back in so the 0s (which are not real 0s) don't influence calculations
		Frs[,,2]= as.matrix(Fr2.rast)
		#Standardize by the mean
		cc.p = cc.p [floor(np.xx/2):floor((np.xx+np.xx/2)-1),floor(np.yy/2):
			floor((np.yy+np.yy/2)-1) ]
		cc.p[is.na(Frs[,,2])]=NA

		nr2_eq_all[ is.na(Frs[,,2])] = NA

		#Equilibrium average competitive fitness
		l1.p[s,rs]=mean(Frs[,,1],na.rm=T)
		y.p[s,rs]=mean(nr2_eq_all,na.rm=T)
		cr_mean.p[s,rs]=mean(cc.p,na.rm=T)
		D.p[s,rs]= cr_mean.p[s,rs]
		
		#Non-linear competitive variance and the spatial storage effect
		#These terms represent the perturbation terms from Snyder 2008 
		var_mu_Us.p[s,rs] = var( c((cc.p)/mean(cc.p,na.rm=T)-1),na.rm=T )
		cov_e_mu_Us.p[s,rs] = cov( c(Frs[,,1]/mean(Frs[,,1],na.rm=T)-1), c( (cc.p)/mean((cc.p),na.rm=T)-1),use="complete.obs" )
		
		#Components of the IGR
		Elam1.p[s,rs]=l1.p[s,rs]*(1/D.p[s,rs]+(cr_mean.p[s,rs]*var_mu_Us.p[s,rs])/D.p[s,rs]^3
			-(cr_mean.p[s,rs]*cov_e_mu_Us.p[s,rs])/D.p[s,rs]^2)+sr[s]-1
		Elam2.p[s,rs]=(l1.p[s,rs]*(1/D.p[s,rs]) +sr[s]-1)^2+2*(l1.p[s,rs]*(1/D.p[s,rs]) +sr[s]-1)*
			(var_mu_Us.p[s,rs]-cov_e_mu_Us.p[s,rs])/D.p[s,rs]^4
		gr1.n.p[s,rs] = exp(Elam1.p[s,rs]-0.5*Elam2.p[s,rs])
		
		#The fitness-density covariance 
		tryCatch( { cov.tmp = get.fd.cov2D(Frs[,,1], cc.p, l1.p[s,rs], D.p[s,rs], gr1.n.p[s,rs], sr[s], a_rr[s])
								cov_lam_vc.p[s,rs]= cov.tmp	[1,1]	
								#Look for anisotropy: that is, look for a difference in cov 
								#between the x and y directions
								sd.tmp = max(cov.tmp[2,1],cov.tmp[2,2]) #Which SD is larger
								if ( abs(cov.tmp[1,2] - cov.tmp[1,1]) >= sd.tmp ) {anis.yes[s,rs]=1}
					} , error=function(e){} )
		
		#The full IGR
		gr1.p[s,rs] = gr1.n.p[s,rs]+cov_lam_vc.p[s,rs]

		save(file="spatial_coexistence_disp_norm2_BCI2.var", "l1.p","y.p", "cr_mean.p","D.p","var_mu_Us.p","cov_e_mu_Us.p",
			"cov_lam_vc.p","Elam1.p", "Elam2.p", "gr1.n.p", "gr1.p")

	}
}

