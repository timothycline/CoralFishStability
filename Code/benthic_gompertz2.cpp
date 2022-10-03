#include <TMB.hpp>
#include <fenv.h>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //DATA
  DATA_MATRIX(obs); // benthic cover observations
  DATA_ARRAY(Covars); //temporal covariate array (1 site x time matrix per covariate)
  DATA_IVECTOR(Habs); //habitat index (must be integer)
  //DATA_VECTOR(Source); //source index
  DATA_INTEGER(K); // Lag for autocorrelation in benthic cover (typically 1)
  DATA_INTEGER(Kbeta); //Lag for covariate effect (testing for longer-range effects of fish community)
  DATA_INTEGER(Kmax); //Maximum lag (need to do lag selection using AICc)
  DATA_INTEGER(Bt); //Are unlagged covariates included (boolean)
  DATA_VECTOR(R_Covar); // Mean fish attribute for recovery rate
  
  //FIXED AND LATENT PARAMETERS
  PARAMETER_VECTOR(R_hab); // Recovery rates by habitat
  PARAMETER(lSigma_R); //Variance of recovery rate regression
  PARAMETER_VECTOR(R_sh); //Recovery rate by site and habitat
  PARAMETER_VECTOR(Beta_Rh); // Effect of fish on recovery rate
  PARAMETER_ARRAY(Beta_hab); //Habitat level effects of fish on benthic cover

  PARAMETER_MATRIX(Ne); //Initial latent benthic cover
  
  
  //PARAMETER_ARRAY(lSigma_beta_hab); //Variance in fish effects across sites
  //PARAMETER_ARRAY(Beta1); //Site level fish effects
  PARAMETER_VECTOR(Phi_h); //Habitat level density dependendence
  PARAMETER_VECTOR(logsdPro); //Process error
  PARAMETER(logsdObs); //ObservationError error
  
  //INDEX VARIABLES
  int timeSteps=obs.row(0).size();
  int nTimeSeries=obs.col(0).size();
  int Ncovar = Covars.cols();
  
  //PARAMETER TRANSFORMATIONS
  Type sdObs=exp(logsdObs);
  vector<Type> sdPro=exp(logsdPro);
  Type Sigma_R = exp(lSigma_R);
  
  Type ans=0;
  
  // Habitat specific density-dependence parameter (deprecated as we fit fixed habitat effects)
  for(int i=0; i<3; i++){
    //ans -= dnorm(Phi_hab(i),Phi,Sigma_Phi,1);
  }
  
  for(int i=0;i<nTimeSeries;i++){
	  
	  //Likelihood for initial values
	  //Initialize values up to lag K
	  for(int l=0;l<Kmax;l++){
	  	//ans -= dnorm(Ne(i,K-(l+1)),Type(0.0),sdObs(i),true);
	  	ans -= dnorm(Ne(i,l),Type(0.0),Type(3.0),1);
	    ans -= dnorm(obs(i,l),Ne(i,l),sdObs,1);
    }
	  
	  //Recovery rate in each site
	  ans -= dnorm(R_sh(i),R_hab(Habs(i)) + Beta_Rh(Habs(i)) * R_Covar(i),Sigma_R,1);
	  
	  Type Pred1;
	  Pred1 = 0.0;
	  //processes occuring at the time step scale
	  for(int n=Kmax; n<timeSteps;n++){	
		  //Ne(i,n) = R_sh(i) + Phi_h(Habs(i)) * exp(Ne(i,n-K)); //DOUBLE TRIPLE CHECK INDEXING HERE
		  Pred1 = R_sh(i) + (1-Phi_h(Habs(i))) * (Ne(i,n-K)); //DOUBLE TRIPLE CHECK INDEXING HERE
		  if(Bt==1){ //if there are unlagged covariates
			for(int B=0;B<Ncovar;B++){
		  		Pred1 += Beta_hab(Habs(i),0,B) * Covars.col(B).col(n).vec()(i); //effect of unlagged covars in 0 position
		  	  	if(Kbeta>0){ // lagged covariate effects
		  			for(int bl=1;bl<(Kbeta+1);bl++){
		  				Pred1 += Beta_hab(Habs(i),bl,B) * Covars.col(B).col(n-bl).vec()(i);
		  			} //bl
		  	  	}
			} //B
		  }else{
			for(int B=0;B<Ncovar;B++){
		  	  	if(Kbeta>0){ ///if there are lagged covariate effects
		  			for(int bl=1;bl<(Kbeta+1);bl++){
		  				Pred1 += Beta_hab(Habs(i),bl-1,B) * Covars.col(B).col(n-bl).vec()(i);
		  			} //bl
		  	  	}
			} //B
		  }
	
		  //regression errors
		  //Ne(i,n) = Pred1;
		  ans -= dnorm(Ne(i,n),Pred1,sdPro(i),1);
		  ans -= dnorm(obs(i,n),Ne(i,n),sdObs,1);
	}	//n	   
  }// i 
  
  REPORT(R_hab); // Recovery rates by habitat
  REPORT(lSigma_R); //Variance of recovery rate regression
  REPORT(R_sh); //Recovery rate by site and habitat
  REPORT(Beta_Rh); // Effect of fish on recovery rate
  REPORT(Ne); //latent benthic cover
  REPORT(Beta_hab); //Habitat level effects of fish on benthic cover
  REPORT(Phi_h); //Habitat level density dependendence
  REPORT(logsdPro); //Observation error
  REPORT(logsdObs); //Observation error
  
  ADREPORT(R_hab); // Recovery rates by habitat
  ADREPORT(lSigma_R); //Variance of recovery rate regression
  ADREPORT(R_sh); //Recovery rate by site and habitat
  ADREPORT(Beta_Rh); // Effect of fish on recovery rate
  ADREPORT(Ne); //latent benthic cover
  ADREPORT(Beta_hab); //Habitat level effects of fish on benthic cover
  ADREPORT(Phi_h); //Habitat level density dependendence
  ADREPORT(logsdPro); //Observation error
  ADREPORT(logsdObs); //Observation error
  
  return ans;
}
