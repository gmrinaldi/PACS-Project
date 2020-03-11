################################
#### TEST SCRIPT ###############
################################

library(fdaPDE)

####### 2D ########

#### square 2D (basic case) ####

rm(list=ls())
graphics.off()

data(square2Ddata)

mesh=create.mesh.2D(nodes=nodes)
# x11()
plot(mesh)
# axis(1)
# axis(2)

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function

set.seed(5847947)

a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
  
}

# Exact solution (pointwise at nodes)
sol_exact=rep(0,dim(mesh$nodes)[1])
for(i in 1: dim(mesh$nodes)[1])
  sol_exact[i]=z(mesh$nodes[i,])

ran=range(sol_exact) 

image(FEM(sol_exact, FEMbasis))

# Set smoothing parameter

lambda= seq(10^-6,10^-3,by=5*10^-6)

GCVFLAG=T
GCVMETHODFLAG='Exact'

data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*(ran[2]-ran[1]))

output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                            GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))


#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####

rm(list=ls())

# Load the mesh
data(simpleMesh2Ddata)

mesh=create.mesh.2D(nodes=nodes, triangles = triangles, order=2)

plot(mesh)

# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh)

set.seed(5847947)

# Exact solution
data_exact=sin(pi*mesh$nodes[,1])

# Plot exact solution
plot(FEM(data_exact,FEMbasis = FEMbasis))

# Locations different from nodes
xobs=runif(min=-0.5,max=0.5,n=80)
yobs=runif(min=-0.5,max=0.5,n=80)
loc=cbind(xobs,yobs)

# Exact data - locations different from nodes
data_exact = sin(pi*xobs)

# Perturbed data - locations different from nodes
data = data_exact + rnorm(n = length(xobs), sd = 0.1)

# Set a vector of smoothing coefficients
lambda = c(10^-4, 1, 10^4)

GCVFLAG=TRUE
GCVMETHODFLAG='Exact'   

# Set PDE parameters
PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)

# deterministic covariate - Nodes locations
cov1_nod=sin(pi*mesh$nodes[,1])

#plot covariate 
image(FEM(cov1_nod,FEMbasis))

# Covariates - Locations different from nodes
cov1_nonod=sin(pi*xobs)
cov2_nonod=rnorm(mean=0, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)

output_CPP = smooth.FEM(observations = data, locations=loc, covariates=W_nonod,
                               FEMbasis = FEMbasis, lambda = lambda, 
                               PDE_parameters = PDE_parameters_anys,
                               GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
#image(output_CPP4$fit.FEM)
output_CPP$beta


#### charotid 2D (space-varying + Boundary Conditions + location not at nodes) ####

rm(list=ls())

data(charotid2Ddata)

plot(mesh)

FEMbasis = create.FEM.basis(mesh)

# Set BC 
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

lambda = 10^-2

set.seed(5839745)

Data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))

GCVFLAG=T
GCVMETHODFLAG='Exact'

Sol = smooth.FEM(locations = SpacePoints,
                       observations = Data, 
                       FEMbasis = FEMbasis, 
                       lambda = lambda, 
                       BC = BC, GCV=GCVFLAG)
image(Sol$fit.FEM)  

plot(Sol$fit.FEM)


#### circular mesh (space-varying + locations at mesh nodes + BC + forcing term) ####

rm(list=ls())

data(peak2Ddata)

# create and plot the triangulation
mesh <- create.mesh.2D(nodes, order = 1)
Tri <- mesh$triangles

plot(mesh,asp=1)

basisobj <- create.FEM.basis(mesh)

# original surface
fdobjtrue <- FEM(data,basisobj)
plot.FEM(fdobjtrue)

# homogeneous dirichlet boundary conditions
BC <- list(BC_indices=which(nodes[,1]^2+nodes[,2]^2>24.5),
           BC_values=rep(0,length(which(nodes[,1]^2+nodes[,2]^2>24.5))))

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

# forcing function
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,1,0)
  output
}
xgrid=seq(from=-5,to=5,by=0.1)
ygrid=seq(from=-5,to=5,by=0.1)
image(xgrid,ygrid,matrix(u_func(expand.grid(xgrid,ygrid)),nrow=length(xgrid),ncol=length(ygrid),byrow=FALSE),main='forcing function',asp=1)

# anisotropy matrix
sigma <- matrix(c(5/6,0,0,1/6),nrow=2,ncol=2,byrow=TRUE)


K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = sigma
  output
}

PDE_parameters <- list(K = K_func, b = b_func, c = c_func, u = u_func)

smoothing <- smooth.FEM(observations=data, FEMbasis=basisobj, lambda= 10^c(-3,-1,1,3,10), PDE_parameters = PDE_parameters, BC = BC)
plot.FEM(smoothing$fit.FEM)

# see the difference with u =0

set.seed(1)

u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
}
PDE_parameters <- list(K = K_func, b = b_func, c = c_func, u = u_func)

smoothing <- smooth.FEM(observations=data, FEMbasis=basisobj, lambda= 10^c(-3,-1,1,3,10), PDE_parameters = PDE_parameters, BC = BC)
plot.FEM(smoothing$fit.FEM)


#### areal square 2D ####

rm(list=ls())

data(square2DarealData)

# Test function

set.seed(5847947)

a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
  
}

# Exact solution (pointwise at nodes)
sol_exact=rep(0,dim(mesh$nodes)[1])
for(i in 1: dim(mesh$nodes)[1])
  sol_exact[i]=z(mesh$nodes[i,])

ran=range(sol_exact)

FEMbasis=create.FEM.basis(mesh)

vals=numeric(dim(mesh$nodes)[1])
valscov=numeric(dim(mesh$nodes)[1])
cov=numeric(dim(mesh$nodes)[1])

for(i in 1:dim(mesh$nodes)[1]){
  vals[i]=obs_areal[RDD[i]]
  valscov[i]=obs_areal[RDD[i]]+1*cov_areal[RDD[i]]
  cov[i]=cov_areal[RDD[i]]
}

image(FEM(vals,FEMbasis))
image(FEM(valscov,FEMbasis))
image(FEM(cov,FEMbasis))
image(FEM(RDD,FEMbasis))

W_areal=cbind(cov_areal)

beta_exact=c(1.2)

ran2=range(obs_areal + W_areal%*%beta_exact)

lambda=10^-6

GCVFLAG=T
GCVMETHODFLAG='Exact'

data = obs_areal + W_areal%*%beta_exact + rnorm(RDD_groups, mean=0, sd=0.05*(ran2[2]-ran2[1]))

smooth_areal<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, covariates = W_areal,
                                incidence_matrix=incidence_matrix, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

#image(smooth_areal2$fit.FEM)

#### FPCA - Kfold cross-validation, locations at nodes ####

rm(list=ls())

data(simpleMesh2Ddata)

mesh<-create.mesh.2D(nodes = nodes,triangles = triangles)

FEMbasis=create.FEM.basis(mesh)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc3=FEM(cos(2*pi*mesh$nodes[,2]), FEMbasis)

plot(eigenfunc1)
plot(eigenfunc2)
plot(eigenfunc3)

truedatarange<-max(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))-min(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))
truecoeff<-cbind(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff)

set.seed(5847947)

nSamples=50

sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.03 
sd_error<-0.05

score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)
score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange)

datamatrix.pointwise.exact<-matrix(score1)%*%t(matrix(eigenfunc1$coeff))+matrix(score2)%*%t(matrix(eigenfunc2$coeff))+matrix(score3)%*%t(matrix(eigenfunc3$coeff))
dm.pointwise.centred.exact<-datamatrix.pointwise.exact-matrix(apply(datamatrix.pointwise.exact,2,mean),ncol=ncol(datamatrix.pointwise.exact),nrow=nrow(datamatrix.pointwise.exact),byrow=TRUE)

lambda=10^c(-6,-5,-4,-3,-2)

validation='KFold'
nfolds=5
GCVMETHOD='Stochastic'
nnodes=dim(mesh$nodes)[1]

error<-rnorm(n=nSamples*nnodes,sd=sd_error*truedatarange)

datamatrix.pointwise<-datamatrix.pointwise.exact+error
dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)

sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred,FEMbasis=FEMbasis,lambda=lambda,nPC=3,validation=validation,GCVmethod=GCVMETHOD, NFolds = nfolds)

plot(sol.pointwise$loadings.FEM)

####### 2.5D #######

rm(list=ls())

#### hub (covariates) ####

data(hub25Ddata)

mesh <- create.mesh.2.5D(nodes = nodes,triangles = triangles)

FEMbasis <- create.FEM.basis(mesh)

set.seed(5847947)

# Locations at nodes
nodesLocations=mesh$nodes

# Exact data - Locations at nodes
nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)
for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*nodesLocations[i+1,1]) +  a2* sin(2*pi*nodesLocations[i+1,2]) +  a3*sin(2*pi*nodesLocations[i+1,3]) +1
}

# Plot the exact solution
plot(FEM(coeff=func_evaluation, FEMbasis = FEMbasis))

lambda=seq(10^-4,10^-3, 0.0001)

cov1=4*sin(2*pi*nodesLocations[,1])*cos(2*pi*nodesLocations[,2])
cov2=rnorm(nnodes, mean=3, sd=0.1)
W=cbind(cov1,cov2)

plot(FEM(coeff = cov1, FEMbasis = FEMbasis))
plot(FEM(coeff = cov2, FEMbasis = FEMbasis))

# Fix betas
beta_exact=c(0.45,0.3) 

ran=range(W%*%beta_exact + func_evaluation)

# Plot exact solution
plot(FEM(W%*%beta_exact + func_evaluation, FEMbasis))

GCVFLAG=T
GCVMETHODFLAG='Exact'

data = func_evaluation + W%*%beta_exact +rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))

output_CPP = smooth.FEM(observations = data, covariates = W,
                               FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG, 
                               GCVmethod = GCVMETHODFLAG)

plot(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

#### areal hub ####

rm(list=ls())

data(hub25DarealData)

set.seed(5847947)

nodesLocations=mesh$nodes

nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*nodesLocations[i+1,1]) +  a2* sin(2*pi*nodesLocations[i+1,2]) +  a3*sin(2*pi*nodesLocations[i+1,3]) +1
}
ran=range(func_evaluation)

FEMbasis=create.FEM.basis(mesh)

plot(FEM(func_evaluation,FEMbasis))

vals=numeric(dim(mesh$nodes)[1])
valscov=numeric(dim(mesh$nodes)[1])
cov=numeric(dim(mesh$nodes)[1])

for(i in 1:dim(mesh$nodes)[1]){
  vals[i]=obs_areal[RDD[i]]
  valscov[i]=obs_areal[RDD[i]]+1*cov_areal[RDD[i]]
  cov[i]=cov_areal[RDD[i]]
}

plot(FEM(vals,FEMbasis))
plot(FEM(valscov,FEMbasis))
plot(FEM(cov,FEMbasis))
plot(FEM(RDD,FEMbasis))

GCVFLAG=TRUE
GCVMETHODFLAG='Exact'

plot(FEM(func_evaluation,FEMbasis))

sol_exact=func_evaluation

W_areal=cbind(cov_areal)
W=cbind(cov1)

beta_exact=c(1)

range(W_areal%*%beta_exact)
ran=range(obs_areal + W_areal%*%beta_exact)

lambda=seq(1.8*10^-5,2*10^-5,by=10^-7)

data = obs_areal + W_areal%*%beta_exact + rnorm(RDD_groups, mean=0, sd=0.05*(ran[2]-ran[1]))

smooth_areal<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, covariates = W_areal,
                                incidence_matrix=incidence_matrix, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

plot(FEM(smooth_areal$fit.FEM$coeff[,which.min(smooth_areal$GCV)],FEMbasis))
plot(FEM(smooth_areal$fit.FEM$coeff[,which.min(smooth_areal$GCV)]+W%*%smooth_areal$beta[,which.min(smooth_areal$GCV)],FEMbasis))

#### FPCA - hub (locations at nodes, stochastic GCV) ####

rm(list=ls())

data(hub25Ddata)

mesh <- fdaPDE::create.mesh.2.5D(nodes = nodes,triangles = triangles)
plot(mesh)

nnodes = mesh$nnodes

FEMbasis=create.FEM.basis(mesh)

set.seed(5847947)

eigenfunc1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)
eigenfunc3=FEM(cos(2*pi*mesh$nodes[,2]), FEMbasis)

plot(eigenfunc1)
plot(eigenfunc2)
plot(eigenfunc3)

truedatarange<-max(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))-min(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))
truecoeff<-cbind(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff)

nSamples=50


sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.025#0.0125
sd_error<-0.05

score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)
score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange)

datamatrix.pointwise.exact<-matrix(score1)%*%t(matrix(eigenfunc1$coeff))+matrix(score2)%*%t(matrix(eigenfunc2$coeff))+matrix(score3)%*%t(matrix(eigenfunc3$coeff))
dm.pointwise.centred.exact<-datamatrix.pointwise.exact-matrix(apply(datamatrix.pointwise.exact,2,mean),ncol=ncol(datamatrix.pointwise.exact),nrow=nrow(datamatrix.pointwise.exact),byrow=TRUE)

# Select validation method
validation="GCV"

# Set lambda
if(is.null(validation)){
  lambda=10^-3
}else if(validation=='GCV' | validation=='KFold'){
  lambda=10^c(-6,-5,-4,-3)
}

# Set GCV options
#GCVMETHOD='Exact'
GCVMETHOD='Stochastic'

error<-rnorm(n=nSamples*nnodes,sd=sd_error*truedatarange)
  
datamatrix.pointwise<-datamatrix.pointwise.exact+error
dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)
  
sol.pointwise<-FPCA.FEM(datamatrix=dm.pointwise.centred, FEMbasis=FEMbasis, lambda=lambda, nPC=3, validation=validation, GCVmethod=GCVMETHOD)
  
plot(sol.pointwise$loadings.FEM)

####### 3D ########


#### sphere 3D (covariates + locations not at nodes + stochastic GCV) ####
  
rm(list=ls())

GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic' #for stochastic GCV (default)
  
  # Build mesh: Sphere
data(sphere3Ddata)
sphere3D<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
plot(sphere3D)
FEMbasis <- create.FEM.basis(sphere3D)
nodesLocations=sphere3D$nodes
  
set.seed(5847947)
  
  # Exact test function
nnodes = sphere3D$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
  
func_evaluation = numeric(nnodes)
  
for (i in 0:(nnodes-1)){
    func_evaluation[i+1] = a1* sin(2*pi*sphere3D$nodes[i+1,1]) +  a2* sin(2*pi*sphere3D$nodes[i+1,2]) +  a3*sin(2*pi*sphere3D$nodes[i+1,3]) +1
  }
ran=range(func_evaluation) 
  
plot(FEM(func_evaluation,FEMbasis))
  # Set smoothing parameter
lambda=c(10^-2)
  
set.seed(5847947)
  
  # Generate locations
nloc = 1000
loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T)
  
ind=NULL
for(row in 1:nloc){
    normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
    if(normvec>0.975)   # check points outside the sphere and remove them
      ind = c(ind,row)  
  }
  
loc=loc[-ind,]
nloc=dim(loc)[1]

  
  # Exact test function - locations different from nodes
func_evaluation2=numeric(nloc)
  for (i in 0:(nloc-1)){
    func_evaluation2[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) +1
  }
ran2=range(func_evaluation2) 
  
set.seed(5847947)
  
cov1_nonod=4*sin(2*pi*loc[,1])+6*sin((2*pi*loc[,2])^2)
cov2_nonod=cos(-2*pi*loc[,3])
  
W2=cbind(cov1_nonod,cov2_nonod)
  
beta_exact=c(0.3,0.5)
  
ran=range(W2%*%beta_exact + func_evaluation2) 
  
data=func_evaluation2 + W2%*%beta_exact + rnorm(nloc,mean=0,sd=0.05*(ran[2]-ran[1]))
  
output_CPP = smooth.FEM(observations = data,locations=loc, covariates=W2,
                                 FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG)
  
plot(output_CPP$fit.FEM)

#### areal sphere 3D ####
  
rm(list=ls())
  
data("sphere3DarealData")
  
FEMbasis=create.FEM.basis(mesh)

set.seed(5847947)
  
# Exact test function
nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
  
func_evaluation = numeric(nnodes)
  
for (i in 0:(nnodes-1)){
    func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
  }
ran=range(func_evaluation)  
  
  
vals=numeric(dim(mesh$nodes)[1])
valscov=numeric(dim(mesh$nodes)[1])
cov=numeric(dim(mesh$nodes)[1])
  
for(i in 1:dim(mesh$nodes)[1]){
    vals[i]=obs_areal[RDD[i]]
    valscov[i]=obs_areal[RDD[i]]+0.8*cov_areal[RDD[i]]
    cov[i]=cov_areal[RDD[i]]
}
  
plot(FEM(vals,FEMbasis))
plot(FEM(valscov,FEMbasis))
plot(FEM(cov,FEMbasis))
plot(FEM(RDD,FEMbasis))

  # Set GCV options
GCVFLAG=TRUE
  #GCVMETHODFLAG='Exact'
GCVMETHODFLAG='Stochastic'
  
  # Set smoothing parameter
lambda=seq(0.000001,0.000009,0.000001)
  
W_areal=cbind(cov_areal)
W=cbind(cov1)
  
beta_exact=c(0.8)
  
plot(FEM(func_evaluation+W%*%beta_exact,FEMbasis))
  
range(W_areal%*%beta_exact) 
ran=range(obs_areal + W_areal%*%beta_exact)
  
data = obs_areal + W_areal%*%beta_exact + rnorm(RDD_groups, mean=0, sd=0.05*(ran[2]-ran[1]))
  
smooth_areal<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, covariates = W_areal,
                                  incidence_matrix=incidence_matrix, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
  
  
#### FPCA - sphere  (locations different from nodes, exact GCV) ####
  
rm(list=ls())

data(sphere3Ddata)
mesh<-create.mesh.3D(sphere3Ddata$nodes,sphere3Ddata$tetrahedrons)
FEMbasis<-create.FEM.basis(mesh)
  
  ## SET PARAMETERS
lambda=10^c(-6,-5,-4,-3)

validation='GCV'
GCVMETHOD='Exact'
  
seed<-1234
nSamples<-50
sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.025
sd_error<-0.05
set.seed(seed)
  
eigenf1=FEM(sin(2*pi*mesh$nodes[,1]), FEMbasis)
eigenf2=FEM(cos(2*pi*mesh$nodes[,1]), FEMbasis)
eigenf3=FEM(cos(2*pi*mesh$nodes[,2]), FEMbasis)
  
plot(eigenf1)
plot(eigenf2)
plot(eigenf3)
  # 
truedatarange.loc<-max(c(eigenf1$coeff,eigenf2$coeff,eigenf3$coeff))-min(c(eigenf1$coeff,eigenf2$coeff,eigenf3$coeff))
truecoeff<-cbind(eigenf1$coeff,eigenf2$coeff,eigenf3$coeff)
  
nloc = 1000
loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T) # randomly generated points
  
ind=NULL
for(row in 1:nloc){
    normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
    if(normvec>0.975)   # remove points that fall outside the sphere
      ind = c(ind,row)
}

loc=loc[-ind,]
nloc=dim(loc)[1]
  
eigenfunc1.loc=sin(2*pi*loc[,1])
eigenfunc2.loc=cos(2*pi*loc[,1])
eigenfunc3.loc=cos(2*pi*loc[,2])
  
score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange.loc)
score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange.loc)
score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange.loc)
  
datamatrix.pointwise.exact.loc<-matrix(score1)%*%t(matrix(eigenfunc1.loc))+matrix(score2)%*%t(matrix(eigenfunc2.loc))+matrix(score3)%*%t(matrix(eigenfunc3.loc))
dm.pointwise.centred.exact.loc<-datamatrix.pointwise.exact.loc-matrix(apply(datamatrix.pointwise.exact.loc,2,mean),ncol=ncol(datamatrix.pointwise.exact.loc),nrow=nrow(datamatrix.pointwise.exact.loc),byrow=TRUE)
  
datamatrix.pointwise.exact2<-matrix(score1)%*%t(matrix(eigenf1$coeff))+matrix(score2)%*%t(matrix(eigenf2$coeff))+matrix(score3)%*%t(matrix(eigenf3$coeff))
dm.pointwise.centred.exact2<-datamatrix.pointwise.exact2-matrix(apply(datamatrix.pointwise.exact2,2,mean),ncol=ncol(datamatrix.pointwise.exact2),nrow=nrow(datamatrix.pointwise.exact2),byrow=TRUE)
error<-rnorm(n=nSamples*dim(loc)[1],sd=sd_error*truedatarange.loc)
  
datamatrix.pointwise.loc<-datamatrix.pointwise.exact.loc+error
dm.pointwise.centred.loc<-datamatrix.pointwise.loc-matrix(apply(datamatrix.pointwise.loc,2,mean),ncol=ncol(datamatrix.pointwise.loc),nrow=nrow(datamatrix.pointwise.loc),byrow=TRUE)
  
sol.pointwise.loc<-FPCA.FEM(locations=loc,datamatrix=dm.pointwise.centred.loc,FEMbasis=FEMbasis,lambda=lambda,nPC=3,validation=validation, GCVmethod=GCVMETHOD)
  
plot(sol.pointwise.loc$loadings.FEM)
  
  
#### example RMSE - hub 2.5D ####  

rm(list=ls())

data(hub25Ddata)

mesh <- create.mesh.2.5D(nodes = nodes,triangles = triangles)

FEMbasis <- create.FEM.basis(mesh)

nodesLocations=mesh$nodes

GCVFLAG=T
GCVMETHODFLAG='Stochastic'

### compute the exact function generating data ###
set.seed(5847947)

nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*nodesLocations[i+1,1]) +  a2* sin(2*pi*nodesLocations[i+1,2]) +  a3*sin(2*pi*nodesLocations[i+1,3]) +1
}
ran=range(func_evaluation)


# Plot exact solution
plot(FEM(coeff=func_evaluation, FEMbasis = FEMbasis))

# Set smoothing parameter
lambda=seq(10^-4,10^-3, 0.0001)


# Set number of simulation trials
N=50 

############# NO COVARIATES, LOCATIONS AT NODES

MSE_nocov=NULL
GCV_param1 = matrix(nrow=N, ncol=3)
colnames(GCV_param1) = c("edf", "stderr", "GCV")
selected_lambda=rep(0,N)
for(i in 1:N){
  
  data = func_evaluation+rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))
  
  output_CPP1 = smooth.FEM(observations = data,
                           FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG, 
                           GCVmethod = GCVMETHODFLAG)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=which.min(output_CPP1$GCV)
  MSE=sum((func_evaluation-output_CPP1$fit.FEM$coeff[,selected_lambda[i]])^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
  
  if(GCVFLAG==TRUE){
    ##STORE GCV PARAMETERS
    GCV_param1[i,1]=output_CPP1$edf[selected_lambda[i]]
    GCV_param1[i,2]=output_CPP1$stderr[selected_lambda[i]]
    GCV_param1[i,3]=output_CPP1$GCV[selected_lambda[i]]
  }
}
RMSE_nocov=sqrt(MSE_nocov)

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)

# Mean RMSE
MSE_hat=mean(MSE_nocov)

plot(FEM(coeff=output_CPP1$fit.FEM$coeff[,selected_lambda[i]], FEMbasis = FEMbasis))


cov1=sin(2*pi*nodesLocations[,1])*cos(2*pi*nodesLocations[,2])
cov2=rnorm(nnodes, mean=0, sd=0.1)
W=cbind(cov1,cov2)

### with covariates 

# Plot the covariates
plot(FEM(coeff = cov1, FEMbasis = FEMbasis))
plot(FEM(coeff = cov2, FEMbasis = FEMbasis))

# Fix betas
beta_exact=c(3,0.5) 

ran2=range(W%*%beta_exact + func_evaluation)

# Plot exact solution
plot(FEM(W%*%beta_exact + func_evaluation, FEMbasis))

# Set smoothing parameter
lambda=seq(10^-4,10^-3, 0.0001)


selected_lambda2=rep(0,N)
MSE_cov_np=NULL  # MSE on non parametric part of the model
MSE_cov_global=NULL  # global MSE
GCV_param2 = matrix(nrow=N, ncol=3)
colnames(GCV_param2) = c("edf", "stderr", "GCV")
betamat=matrix(data=NA, nrow = 2, ncol=N)

for(i in 1:N){
  
  data2 = func_evaluation + W%*%beta_exact +rnorm(nnodes,mean=0,sd=0.05*(ran2[2]-ran2[1]))
  
  output_CPP2 = smooth.FEM(observations = data2, covariates = W,
                           FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG, 
                           GCVmethod = GCVMETHODFLAG)
  
  selected_lambda2[i]=which.min(output_CPP2$GCV)
  
  #plot(output_CPP2$fit.FEM)
  #plot(FEM(output_CPP2$fit.FEM$coeff+W%*%output_CPP2$beta,FEMbasis))
  
  MSE_np=sum((func_evaluation-output_CPP2$fit.FEM$coeff[,selected_lambda2[i]])^2)/nnodes
  MSE_global=sum((func_evaluation+W%*%beta_exact-output_CPP2$fit.FEM$coeff[,selected_lambda2[i]]-W%*%output_CPP2$beta[,selected_lambda2[i]])^2)/nnodes
  
  MSE_cov_np=c(MSE_cov_np, MSE_np)
  MSE_cov_global=c(MSE_cov_global, MSE_global)
  
  betamat[,i]=output_CPP2$beta[,selected_lambda2[i]]
  
  if(GCVFLAG==TRUE){
    ##STORE GCV PARAMETERS
    GCV_param2[i,1]=output_CPP2$edf[selected_lambda2[i]]
    GCV_param2[i,2]=output_CPP2$stderr[selected_lambda2[i]]
    GCV_param2[i,3]=output_CPP2$GCV[selected_lambda2[i]]
  }
  
}
# table(selected_lambda2)


RMSE_cov_np=sqrt(MSE_cov_np)
RMSE_cov_global=sqrt(MSE_cov_global)

# Boxplot of RMSE
x11()
boxplot(RMSE_cov_np,RMSE_cov_global, names = c('RMSE_np','RMSE_tot'))


# Mean RMSE
MSE_hat_cov_np=mean(MSE_cov_np) 
sqrt(MSE_hat_cov_np)

MSE_hat_cov_global=mean(MSE_cov_global) 
sqrt(MSE_hat_cov_global)


# Mean std error (computed if GCV=T)
mean(GCV_param2[,2])


plot(FEM(coeff=output_CPP2$fit.FEM$coeff[,selected_lambda2[i]]+W%*%output_CPP2$beta[,selected_lambda2[i]], FEMbasis = FEMbasis))

# 0.05*(ran2[2]-ran2[1])

# Check the quality of the estimation for betas
x11()
par(mfrow=c(1,2))
boxplot(betamat[1,], main='Beta 1')
abline(h=beta_exact[1], col='red')
boxplot(betamat[2,], main='Beta 2')
abline(h=beta_exact[2], col='red')   

