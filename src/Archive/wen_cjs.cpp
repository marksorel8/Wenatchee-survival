#include <TMB.hpp>



template<class Type>
Type objective_function<Type>::operator() ()
{
//~~~~~~~~~~~~~~~~~~~
// Data
//~~~~~~~~~~~~~~~~~~~
  DATA_INTEGER(n_occ); //number of survival/ reacputure occasions;
  DATA_INTEGER(N_unique_CH); //number of survival/ reacputure occasions;
  DATA_MATRIX(CH);     //capture histories (excluding occasion at marking (which we are conditioning on))
  DATA_IVECTOR(freq); //frequency of capture histories
  DATA_MATRIX(X_phi); //fixed effect design matrix for phi
  DATA_MATRIX(X_p);   // fixed effect design matrix for psi
  DATA_IMATRIX(Phi_pim); //index of phi parameter vector for a given Ch x occasion
  DATA_IMATRIX(p_pim);   //index of p parameter vector for a given Ch x occasion  
  
  //~~~~~~~~~~~~~~~~~~~
  // Parameters
  //~~~~~~~~~~~~~~~~~~~
  PARAMETER_VECTOR(beta_phi);  //Phi fixed effect coefficients
  PARAMETER_VECTOR(beta_p);    //p fixed effect coefficients
  
  
  // Joint negative log-likelihood
  parallel_accumulator<Type> jnll(this);
  
  // Linear predictors
  vector<Type> eta_phi = X_phi*beta_phi;
  vector<Type> eta_p = X_p*beta_p;
  
  
  // Apply link
  vector<Type> phi=invlogit(eta_phi);
  vector<Type> p=invlogit(eta_p); 
  REPORT(phi);
  REPORT(p);
  
  
  // Observation likelihood
  
  Type p_alive = 0; //prob being alive
  Type p_dead=0;    //prob being dead
  Type u = 0;
  Type NLL_it=0;
  
  for(int n=0; n<N_unique_CH; n++){ // loop over individual unique capture histories
    p_alive=Type(1); //condition on known alive at capture
    p_dead=Type(0); //known not dead on capture
    NLL_it=Type(0); //initiale capture history NLL at 0
  for(int t=0; t<n_occ; t++){       //loop over occasions (excluding capture occasion)
    
    p_dead += Type((Type(1)-phi(Phi_pim(n,t)))*p_alive); //prob die or stay dead
    p_alive *= Type(phi(Phi_pim(n,t))); //prob stay alive
   
    
    
    p_alive *= Type(p(p_pim(n,t))*CH(n,t)+ (Type(1)-p(p_pim(n,t)))*(Type(1)-CH(n,t))); //prob observation given alive
    p_dead *= Type(Type(1)-CH(n,t)); //prob observation given dead
    
    
    u = p_alive+p_dead;  //sum of probs
    p_alive = p_alive/u; //normalize probs
    p_dead = p_dead/u;   //normalize probs
    NLL_it  +=log(u);    //accumulate nll
  }
  jnll-=(NLL_it*freq(n));  //multiply times frequency of cpature history and subtract from jnll
  }
  
  return(jnll);
}