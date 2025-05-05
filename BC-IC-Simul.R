# R code for data generation

#Case 1: alpha=0

LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}
data_0_BC=function(n,b0,b1,b2,g1,g2,lam,cenrate){ 
  # lam is baseline exponential parameter 
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  x1=rbinom(n=n,size=1,prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  phi=exp(b0+(b1*x1)+(b2*x2))
  U=runif(n,min=0,max=1)
  C=rexp(n,rate=cenrate)
  p0 = exp(-phi)
  count.obs=0
  count.cure=0
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i]=C[i]
      R[i]=Inf
      d[i]=0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y=-exp(-(g1*x1[i])-(g2*x2[i]))
      *(1/lam)*log(1+(exp(-b0-(b1*x1[i])-(b2*x2[i]))
                      *log(p0[i]+((1-p0[i])*U1))))
      t=min(y,C[i])
      if(t==C[i]){    
        L[i]=C[i]
        R[i]=Inf
        d[i]=0
      }else{
        len=runif(1,0.2,0.7)
        l=runif(1,0,1)
        ans=LR_int(t,len,l)
        L[i]=ans[1]
        R[i]=ans[2]
        d[i]=1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  print("count.cure")
  print(count.cure)
  print("count.obs")
  print(count.obs)
  return(data.frame(L,R,d,x1,x2))
}



# Case 2: alpha in (0,1]

LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}
data_gen_BC=function(n,alpha,b0,b1,b2,g1,g2,lam,cenrate){ 
  # alpha is the BC index parameter, lam is baseline exponential parameter 
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  x1=rbinom(n=n,size=1,prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  phi=exp(b0+(b1*x1)+(b2*x2))/(1+(alpha*exp(b0+(b1*x1)+(b2*x2))))
  U=runif(n,min=0,max=1)
  C=rexp(n,rate=cenrate)
  p0 = (1-(alpha*phi))^(1/alpha)
  count.obs=0
  count.cure=0
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i]=C[i]
      R[i]=Inf
      d[i]=0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y=-exp(-(g1*x1[i])-(g2*x2[i]))
      *(1/lam)
      *log(((alpha*phi[i])+((p0[i]+((1-p0[i])*U1))^alpha)-1)/(alpha*phi[i]))
      t=min(y,C[i])
      if(t==C[i]){    
        L[i]=C[i]
        R[i]=Inf
        d[i]=0
      }else{
        len=runif(1,0.2,0.7)
        l=runif(1,0,1)
        ans=LR_int(t,len,l)
        L[i]=ans[1]
        R[i]=ans[2]
        d[i]=1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  print("count.cure")
  print(count.cure)
  print("count.obs")
  print(count.obs)
  return(data.frame(L,R,d,x1,x2))
}


# R code for estimation of parameters (by simultaneous maximization) and SEs

data2=data1 
nt1=nrow(data2)

### Separating the censored data from the uncensored data
inc=which(data2$d==1) 
ic=which(data2$d==0)
ninc=length(inc)
nic=length(ic)
pinc=ninc/nt1
pic=nic/nt1

### Sorting time points excluding the 
#"Inf" time points since this will be required 
#for the choice of cut points
timepts=sort(unique(c(data2$L,data2$R)))
timn=length(timepts)
timepts=timepts[-timn]
plot(timepts,type="l")
hist(timepts)
data3=data2

### Some more pre-processing 
data3$Timept1=data3$L
data3$Timept2=data3$R
data3=data3[,c(4,5,6,7)]

### Defining the X and Z vectors
Xvec=as.matrix(data3[,c(1,2)]) 
Zvec=cbind(1,Xvec)

###################################################
########### All necessary functions ###############
###################################################
#### Some intermediate functions to facilitate calculations 
#### and choice of cut points and hazard values at those points 
ctppsi_func=function(num_lines){
  ctp_choice=c(0,quantile(timepts,seq(0,0.9,0.05)))
  NCP=length(ctp_choice)
  psi_choice=c(0.001,rep(0.1,(NCP-1)))
  for(uu in 2:NCP){
    psi_choice[uu]=(psi_choice[uu-1])+(psi_choice[uu])+(psi_choice[uu])^2
  }
  if(num_lines==1)
  {ctp_cons=ctp_choice[c(1,NCP)];psi_cons=psi_choice[c(1,NCP)]}
  if(num_lines==2)
  {ctp_cons=ctp_choice[c(1,floor(NCP/2),NCP)];
  psi_cons=psi_choice[c(1,floor(NCP/2),NCP)]}
  if(num_lines==3)
  {ctp_cons=ctp_choice[c(1,floor(NCP/3),2*floor(NCP/3),NCP)];
  psi_cons=psi_choice[c(1,floor(NCP/3),2*floor(NCP/3),NCP)]}
  if(num_lines==4)
  {ctp_cons=ctp_choice[c(1,floor(NCP/4),2*floor(NCP/4),
                         3*floor(NCP/4),NCP)];
  psi_cons=psi_choice[c(1,floor(NCP/4),2*floor(NCP/4),3*floor(NCP/4),NCP)]}    
  return(c(ctp_cons,psi_cons))
}    

########################################################################
##########  Defining optimization functions and calculation of std errors
pla_optim=function(num_lines){
  ### Cut points and initial estimates of hazard
  chopts=ctppsi_func(num_lines)[1:(num_lines+1)]
  psinit=ctppsi_func(num_lines)[(num_lines+2):(2*num_lines+2)]
  thetain=round(c(psinit,gaminit,betinit),3)
  ### PLA function
  pla=function(tt,psi,ctp){
    NN=length(psi) 
    hj=rep(0,NN)
    cj=rep(0,NN)
    fj=rep(0,NN)
    sj=rep(0,NN)
    ej=rep(0,NN)
    II=rep(0,NN)
    for(ii in 2:NN){
      sj[ii]=(psi[ii]-psi[ii-1])/(ctp[ii]-ctp[ii-1])
      cj[ii]=psi[ii]-ctp[ii]*sj[ii]
      hj[ii]=cj[ii]+sj[ii]*tt
      ej[ii]=cj[ii]*(min(tt,ctp[ii])-ctp[ii-1])
      fj[ii]=ej[ii]+(sj[ii]/2)*((min(tt,ctp[ii]))^2-(ctp[ii-1])^2)
    }
    for(ii in 2:NN){
      if(tt>=ctp[ii-1]&tt<ctp[ii]){II[ii]=1}
    }
    iii=which(II==1)
    
    if(is.integer(iii) && length(iii) == 0){
      hval=psi[NN]
      Hval=sum(fj[1:NN])
      if(Hval<0){Hval=0}
      Sval=exp(-Hval)
      Fval=1-Sval
      fval=hval*Sval
      return(c(hval,Hval,Sval,Fval,fval))
    }
    else{    
      hval=hj[iii]
      Hval=sum(fj[1:iii])
      if(Hval<0){Hval=0}
      Sval=exp(-Hval)
      Fval=1-Sval
      fval=hval*Sval
      return(c(hval,Hval,Sval,Fval,fval))
    }
  }
  
  ### Complete data log-likelihood function with baseline 
  ### hazard being estimated by PLA 
  comloglik=function(theta){
    alp=alp_0
    psidf=theta[1:(length(psinit))];
    gamth=theta[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta[(length(psinit)+length(gaminit)+1):(length(theta))]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
      valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
    }
    #### alpha=0
    if(alp==0){
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=exp(-thfunvec[ii]*intml[ii])
        Spvecr[ii]=exp(-thfunvec[ii]*intmr[ii])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      } 
    }
    finval=sum(log((Spvecl[inc]-Spvecr[inc])))
    +sum((1-w0)*log((p0vec[ic])))+sum(w0*log(((1-p0vec[ic])*Suvecl[ic])))
    return(-finval)
  }
  
  ### Observed data log-likelihood function with baseline hazard being 
  ### estimated by PLA 
  loglikf=function(theta){
    alp=alp_0
    psidf=theta[1:(length(psinit))];
    gamth=theta[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta[(length(psinit)+length(gaminit)+1):(length(theta))]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
      valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
    }
    #### alpha=0
    if(alp==0){
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=exp(-thfunvec[ii]*intml[ii])
        Spvecr[ii]=exp(-thfunvec[ii]*intmr[ii])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    finval=sum(log(((Spvecl[inc]-Spvecr[inc]))))+sum(log((Spvecl[ic])))
    return(finval)
  }
  
  ### Main Loop and optimization #########
  ########################################
  ctt=0;eps0=0.001;eps=10;theta0=thetain;
  
  while(eps>eps0){
    
    alp=theta0[length(theta0)]
    psidf=theta0[1:(length(psinit))];
    gamth=theta0[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta0[(length(psinit)+length(gaminit)+1):(length(theta0)-1)]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
    
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=exp_f(psidf,data3$Timept1[ii])
      valmatr[ii,]=exp_f(psidf,data3$Timept2[ii])
    }
    
    #### alpha=0
    if(alp==0){
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=exp(-thfunvec[ii]*intml[ii])
        Spvecr[ii]=exp(-thfunvec[ii]*intmr[ii])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    
    w0=((1-p0vec[ic])*Suvecl[ic])/Spvecl[ic]
    
    ######### Optimizing objective function
    #######################################
    obj1=optim(theta0,comloglik_exp,method="Nelder-Mead")
    theta1=obj1$par
    #eps=max(abs((theta1-theta0)/theta0))
    eps=sqrt(sum((theta1-theta0)^2))
    theta0=theta1;
    print(round(theta0,4))
    ctt=ctt+1;
    LLO=loglikf_exp(theta1)
    print(c(LLO,ctt,eps))
  }
  
  estm=theta0
  estm2=round(estm,4)
  
  #### Calculation of SE ######
  #############################
  H1=-(hessian(loglikf_exp,estm))
  detH1=det(H1)
  H2=(pinv(H1))
  hdiag=diag(H2)
  seval=round(sqrt(hdiag),4)
  lcl=estm2-1.96*seval
  ucl=estm2+1.96*seval
  retval=list()
  retval[[1]]=round(estm2,4)
  retval[[2]]=round(seval,4)
  retval[[3]]=round(lcl,4)
  retval[[4]]=round(ucl,4)
  retval[[5]]=round(loglikf_exp(estm),4)
  return(retval)
  
}


