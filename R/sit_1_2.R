#functions to simulate simple MSF and calculate SIT methods 1 and 2
# on Jim scott's spreadsheets
#Author: Catarina Wor
#date: April 2020

#' create_example_input 
#' 
#' 
#'
#' @export
#' @return A list containing the simulation inputs of jim's excel sheet
#' 
#'
create_example_input <- function() {

    ex <- list()

    ex$nfisheries <- 10
    ex$nage <- 3

    ex$msf_flag <- matrix(1,ncol=ex$nage, nrow=ex$nfisheries)
    ex$msf_flag[2,]<-0

    # msf_ release mortality rate
    ex$msfF <- matrix(
	c(0.2000, 0.2000, 0.2000,
      0.1400, 0.1400, 0.1400,
      0.4670, 0.3070, 0.3040,
      0.1700, 0.1290, 0.1350,
      0.1700, 0.1290, 0.1350,
      0.6000, 0.6000, 0.6000,
      0.1700, 0.1290, 0.1350,
      0.1700, 0.1290, 0.1350,
      0.1700, 0.1290, 0.1350,
      0.1700, 0.1290, 0.1350),
	ncol=ex$nage, nrow=ex$nfisheries, byrow=T)

    #Proportion to Ocean Subpopulations
    ex$g <- c(.8, .2)

    #Natural Mortality Rate
    ex$M <- c(.4,.3,.2)

    #Initial Cohort size
    ex$Nu <- 50000
    ex$Nm <- 50000

    #maturation rates
    ex$pa <- c(0.0924,	0.9221,	1.0000)

    #Harvest Rates (Base)
    ex$hr <- matrix(
    c(0.2000,	0.2000,	0.2000,
        0.1000,	0.1000,	0.1000,
        0.0203,	0.0295,	0.0297,
        0.0547,	0.0769,	0.0611,
        0.0059,	0.0086,	0.0066,
        0.0769,	0.0897,	0.0676,
        0.0021,	0.0017,	0.0009,
        0.0007,	0.0010,	0.0007,
        0.0033,	0.0006,	0.0036,
        0.4000,	0.0047,	0.0064),
    ncol=ex$nage, nrow=ex$nfisheries, byrow=T)


    #Discard Mortality Scalar
    ex$dsc <- matrix(
    c(1.9091,	0.1395,	0.1000,
    1.4318,	0.1047,	0.0750,
    0.0500,	0.0705,	0.0163,
    0.4583,	0.3198,	0.1102,
    0.4583,	0.3198,	0.1102,
    0.1121,	0.0789,	0.0159,
    0.2941,	0.3253,	0.1102,
    0.4583,	0.3198,	0.1102,
    0.4583,	0.3198,	0.1102,
    1.0000,	0.3153,	0.0909),
    ncol=ex$nage, nrow=ex$nfisheries, byrow=T)

    return(ex)

}





#' Simulate data with mark selective fisheries 
#' 
#' 
#' @param exdata A list containing data inputs, similar to that
#'  produced as an output of create_example_input
#' @export
#' @return A list containing the simulation outputs of a simple msf fishery
#' 
#'
simpop <- function(exdata) {

    with(exdata,{

    #Cohort before fisheries
    Anot <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    A <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TM <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TC <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TD <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    Tmat <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))


    #Harvest Rates (Base)
    hru <- matrix(NA, ncol=nage, nrow=nfisheries, byrow=T)
    lambda <- matrix(NA, ncol=nage, nrow=nfisheries+1, byrow=T)
    escapement <- matrix(NA, ncol=nage, nrow=2)
    hrm <- hr

    #loop over ages
    for(a in seq_len(nage)){
        #loop over 2 groups marked and unmarked
        for(n in 1:2){
        
            #Ocean fisheries - ocean fisheries happen simulataneously but different proportions of the population are exposed to different set of  oceanfisheries
            for(i in seq_along(g)){

                if(a==1){
                    Anot[[n]][i,a] <- Nm*g[i]         
                }else{
                    Anot[[n]][i,a] <- A[[n]][i,a-1]-
                                     TC[[n]][i,a-1] - 
                                     TD[[n]][i,a-1] - 
                                     Tmat[[n]][i,a-1]
                }
        
                TM[[n]][i,a] <- Anot[[n]][i,a]* M[a]
                A[[n]][i,a] <- Anot[[n]][i,a] - TM[[n]][i,a] 
            
                if(n==1){
                   TC[[n]][i,a] <- A[[n]][i,a] * hr[i,a]
                   TD[[n]][i,a] <- TC[[n]][i,a] * dsc[i,a]
                }else if(n==2){
                    TC[[n]][i,a] <- ifelse(msf_flag[i,a],
                                 A[[n]][i,a] *hr[i,a]*msfF[i,a],
                                 A[[n]][i,a] *hr[i,a])
                    TD[[n]][i,a] <- A[[n]][i,a]*hr[i,a]*dsc[i,a]
                    
                    lambda[i,a] <- A[[2]][i,a]/A[[1]][i,a]

                }
               
                Tmat[[n]][i,a] <- (A[[n]][i,a] - TC[[n]][i,a] - TD[[n]][i,a])* pa[a]
            }

            #Terminal Fisheries
            for(j in seq(length(g)+1, nfisheries)){

                if(j==length(g)+1){
                   A[[n]][j,a] <- sum(Tmat[[n]][1:(j-1),a])
                
                }else{

                    A[[n]][j,a] <- A[[n]][j-1,a]-
                              TC[[n]][j-1,a]-
                              TD[[n]][j-1,a]
                }
        
                if(n==1){
                    TC[[n]][j,a] <- A[[n]][j,a] * hr[j,a]
                    TD[[n]][j,a] <- TC[[n]][j,a] * dsc[j,a]
                }else if(n==2){
                    TC[[n]][j,a] <- ifelse(msf_flag[j,a],
                                 A[[n]][j,a] * hr[j,a] * msfF[j,a],
                                 A[[n]][j,a] * hr[j,a])
                    TD[[n]][j,a] <- A[[n]][j,a]*hr[j,a]*dsc[j,a]
                    hru[j,a] <-  TC[[n]][j,a]/A[[n]][j,a]
                    lambda[j,a] <- A[[2]][j,a]/A[[1]][j,a]  
                }
            }
            escapement[n,a] <- A[[n]][j,a] - TC[[n]][j,a] - TD[[n]][j,a]       
        }
    
        lambda[j+1,a] <- escapement[2,a]/escapement[1,a]
    }
    hrm[1,] <- TC[[1]][1,]/(A[[1]][1,]+A[[1]][2,])
    hrm[2,] <- TC[[1]][2,]/(A[[1]][1,]+A[[1]][2,])
    hru[1,] <- TC[[2]][1,]/(A[[2]][1,]+A[[2]][2,])
    hru[2,] <- TC[[2]][2,]/(A[[2]][1,]+A[[2]][2,])
    
    a <- append(exdata, list(lambda=lambda,
        hru=hru,
        hrm=hrm,
        escapement = escapement,
        Anot = Anot, 
        A = A,
        TM = TM,
        TC = TC, 
        TD = TD,
        Tmat = Tmat))

    return(a)
    })
}






#' Simulate and estimate parameters using sit method 1 
#' 
#' 
#' @param simdata A list containing data inputs, similar to that
#'  produced as an output of create_example_input
#' @export
#' @return A list containing the simulation outputs of a simple msf fishery
#' 
#'
sit1 <- function(simdata) {

	with(simdata,{

    # Estimation

    Nnote <- matrix(NA, ncol=nage, nrow=2)
    Ne <- matrix(NA, ncol=nage, nrow=2)
    TMe <- matrix(NA, ncol=nage, nrow=2)
    Ae <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
    	unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TCe <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
    	unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TDe <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
    	unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    Tmate <- matrix(NA, ncol=nage, nrow=2)
    hre <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
    	unmarked=matrix(NA,nrow=nfisheries,ncol=nage))


    #Harvest Rates (Base)
    #hrue <- matrix(NA, ncol=nage, nrow=nfisheries, byrow=T)
    lambdae<- matrix(NA, ncol=nage, nrow=nfisheries+1, byrow=T)
    escapemente <- matrix(NA, ncol=nage, nrow=2)

    for(a in seq(nage,1)){
        #loop over 2 groups marked and unmarked
        lambdae[nfisheries+1,a] <- lambda[nfisheries+1,a]
        for(n in 1:2){
    
        	if(n==1){
        		escapemente[n,a] <- escapement[n,a]
        	}else{
        		escapemente[n,a] <- lambdae[nfisheries+1,a] * escapemente[n-1,a]
        	}   	
        	#terminal fisheries
        	for(j in seq(nfisheries,1)){
    
        		#key assumption of SIT1 - lambda only gete updates as escapement
        		lambdae[j,a] <- lambdae[j+1,a]
     
                if(n==1){
                	TCe[[n]][j,a] <- TC[[n]][j,a]
                	TDe[[n]][j,a] <- TCe[[n]][j,a] * dsc[j,a]
                }else{
                	TCe[[n]][j,a] <- ifelse(msf_flag[j,a],
                	                 lambdae[nfisheries+1,a]*msfF[j,a]*TC[[n-1]][j,a],
                	                 lambdae[nfisheries+1,a]*TCe[[n-1]][j,a])
                	TDe[[n]][j,a] <- ifelse(msf_flag[j,a],
                	                 dsc[j,a]*TCe[[n]][j,a]/msfF[j,a] ,
                	                 dsc[j,a]*TCe[[n]][j,a])
                }
        		
        		if(j>length(g)){
        		    if(j==nfisheries){ 
        			    Ae[[n]][j,a] <- TCe[[n]][j,a] + TDe[[n]][j,a] + escapemente[n,a]
        		    }else{
        			    Ae[[n]][j,a] <- TCe[[n]][j,a] + TDe[[n]][j,a] + Ae[[n]][j+1,a]
        		    }
        		
        			hre[[n]][j,a] <- TCe[[n]][j,a]/Ae[[n]][j,a]
        		}
        	}
    
        	Tmate[n,a] <- Ae[[n]][length(g)+1,a]
       	
        	if(a==nage){
        		Ne[n,a] <- Ae[[n]][length(g)+1,a] + 
        	    sum(TCe[[n]][1:length(g),a]) + 
        	    sum(TDe[[n]][1:length(g),a])
        	}else{
        		Ne[n,a] <- Ae[[n]][length(g)+1,a] + 
        	    sum(TCe[[n]][1:length(g),a]) + 
        	    sum(TDe[[n]][1:length(g),a]) +
        	    Nnote[n,a+1]
        	}
        	
        	TMe[n,a] <- M[a]/(1-M[a]) * Ne[n,a]
        	Nnote[n,a] <- Ne[n,a] + TMe[n,a]
    
            #note this difference with simulated data
        	hre[[n]][1:length(g),a] <- TCe[[n]][1:length(g),a]/Ne[n,a]
        }
    }


    a <- list(data = simdata,
    simulated = list(Anot = Anot, 
        A = A,
        TM = TM,
        TC = TC, 
        TD = TD,
        Tmat = Tmat),
    estimated = list( Nnote = Nnote,
        Ne = Ne,
        TMe = TMe,
        Ae = Ae,   
        TCe = TCe,
        TDe = TDe,
        Tmate = Tmate, 
        hre = hre)
        )

    return(a)
    })
}



#' Simulate and estimate parameters using sit method 2 
#' 
#' 
#' @param simdata A list containing data inputs, similar to that
#'  produced as an output of create_example_input
#' @export
#' @return A list containing the simulation outputs of a simple msf fishery
#' 
#'
sit2 <- function(simdata) {

    with(simdata,{

    # Estimation

    Nnote <- matrix(NA, ncol=nage, nrow=2)
    Ne <- matrix(NA, ncol=nage, nrow=2)
    TMe <- matrix(NA, ncol=nage, nrow=2)
    Ae <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TCe <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    TDe <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))
    Tmate <- matrix(NA, ncol=nage, nrow=2)
    hre <- list(marked=matrix(NA,nrow=nfisheries,ncol=nage),
        unmarked=matrix(NA,nrow=nfisheries,ncol=nage))


    #Harvest Rates (Base)
    #hrue <- matrix(NA, ncol=nage, nrow=nfisheries, byrow=T)
    lambdae<- matrix(NA, ncol=nage, nrow=nfisheries+1, byrow=T)
    escapemente <- matrix(NA, ncol=nage, nrow=2)

    for(a in seq(nage,1)){
        #
        lambdae[nfisheries+1,a] <- lambda[nfisheries+1,a]
        #loop over 2 groups marked and unmarked
        for(n in 1:2){
    
            if(n==1){
                escapemente[n,a] <- escapement[n,a]
            }else{
                escapemente[n,a] <- lambdae[nfisheries+1,a] * escapemente[n-1,a]
            }       
            #terminal fisheries
            for(j in seq(nfisheries,length(g)+1)){
    
                if(n==1){
                    TCe[[n]][j,a] <- TC[[n]][j,a]
                    TDe[[n]][j,a] <- TCe[[n]][j,a] * dsc[j,a]

                    if(j>length(g)){ #if not terminal
                        if(j==nfisheries){ 
                            Ae[[n]][j,a] <- TCe[[n]][j,a] + TDe[[n]][j,a] + escapemente[n,a]
                        }else{
                            Ae[[n]][j,a] <- TCe[[n]][j,a] + TDe[[n]][j,a] + Ae[[n]][j+1,a]
                        }
                        hre[[n]][j,a] <- TCe[[n]][j,a]/Ae[[n]][j,a]
                    }

                }else{
                    
                    if(j==nfisheries){ 
                        Ae[[n]][j,a] <- ifelse(msf_flag[j,a],
                                     lambdae[nfisheries+1,a]*(Ae[[n-1]][j,a]-TCe[[n-1]][j,a]-TDe[[n-1]][j,a])/
                                     (1-msfF[j,a]*hre[[n-1]][j,a]-dsc[j,a]*hre[[n-1]][j,a]),
                                     lambdae[nfisheries+1,a]*(Ae[[n-1]][j,a]-TCe[[n-1]][j,a]-TDe[[n-1]][j,a])/
                                     (1-msfF[j,a]*hre[[n-1]][j,a]-hre[[n-1]][j,a]))
                    }else{
                        Ae[[n]][j,a] <- ifelse(msf_flag[j,a],
                                     lambdae[j+1,a]*(Ae[[n-1]][j,a]-TCe[[n-1]][j,a]-TDe[[n-1]][j,a])/
                                     (1-msfF[j,a]*hre[[n-1]][j,a]-dsc[j,a]*hre[[n-1]][j,a]),
                                     lambdae[j+1,a]*(Ae[[n-1]][j,a]-TCe[[n-1]][j,a]-TDe[[n-1]][j,a])/
                                     (1-msfF[j,a]*hre[[n-1]][j,a]-hre[[n-1]][j,a]))
                    }

                    TCe[[n]][j,a] <- ifelse(msf_flag[j,a],
                                     Ae[[n]][j,a]*msfF[j,a]*hre[[n-1]][j,a],
                                     lambdae[nfisheries+1,a]*TCe[[n-1]][j,a])
                    TDe[[n]][j,a] <- ifelse(msf_flag[j,a],
                                     dsc[j,a]*TCe[[n]][j,a]/msfF[j,a] ,
                                     dsc[j,a]*TCe[[n]][j,a])
                    lambdae[j,a] <- Ae[[n]][j,a]/Ae[[n-1]][j,a]
                }  
            }

            #ocean Fisheries
            lambdae[1:length(g),a] <- lambdae[length(g) + 1,a]
            if(n==1){
                TCe[[n]][1:length(g),a] <- TC[[n]][1:length(g),a]
                TDe[[n]][1:length(g),a] <- TCe[[n]][1:length(g),a]*dsc[j,a]
            }else{
                TCe[[n]][1:length(g),a] <- ifelse(msf_flag[1:length(g),a],
                    lambdae[1:length(g),a]*msfF[1:length(g),a]*TCe[[n]][1:length(g),a],
                    lambdae[1:length(g),a]*TCe[[n]][1:length(g),a])
                TDe[[n]][1:length(g),a] <- ifelse(msf_flag[1:length(g),a],
                    TCe[[n]][1:length(g),a]/msfF[1:length(g),a]*dsc[1:length(g),a],
                    TCe[[n]][1:length(g),a]*dsc[1:length(g),a])
            }
    
            #parei aqui
            if(a==nage){
                Ne[n,a] <- Ae[[n]][length(g)+1,a] + 
                sum(TCe[[n]][1:length(g),a]) + 
                sum(TDe[[n]][1:length(g),a])
            }else{
                Ne[n,a] <- Ae[[n]][length(g)+1,a] + 
                sum(TCe[[n]][1:length(g),a]) + 
                sum(TDe[[n]][1:length(g),a]) +
                Nnote[n,a+1]
            }
            
            TMe[n,a] <- M[a]/(1-M[a]) * Ne[n,a]
            Nnote[n,a] <- Ne[n,a] + TMe[n,a]

            Tmate[n,a] <- Ae[[n]][length(g)+1,a]/(Ne[n,a] -  (sum(TCe[[n]][1:length(g),a]) + 
                sum(TDe[[n]][1:length(g),a])))
    
            #note this difference with simulated data
            hre[[n]][1:length(g),a] <- TCe[[n]][1:length(g),a]/Ne[n,a]
        }
    }


    a <- list(data = simdata,
    simulated = list(Anot = Anot, 
        A = A,
        TM = TM,
        TC = TC, 
        TD = TD,
        Tmat = Tmat),
    estimated = list( Nnote = Nnote,
        Ne = Ne,
        TMe = TMe,
        Ae = Ae,   
        TCe = TCe,
        TDe = TDe,
        Tmate = Tmate, 
        hre = hre)
        )
    return(a)
    })
}


