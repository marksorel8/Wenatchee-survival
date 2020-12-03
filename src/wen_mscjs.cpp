#include <TMB.hpp>



template<class Type>
struct pim: vector<matrix<int> > {

  pim(SEXP x){ // Constructor
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP m = VECTOR_ELT(x, i);
    (*this)(i) = asMatrix<int>(m);
    }

  }
};



template<class Type>
Type objective_function<Type>::operator() ()
{
//~~~~~~~~~~~~~~~~~~~
// Data
//~~~~~~~~~~~~~~~~~~~
DATA_INTEGER(n_OCC);      //number of total survival/ reacputure occasions;
DATA_INTEGER(nDS_OCC);      //number of downstream survival/ reacputure
DATA_INTEGER(n_states);      //number of possible adult return ages (i.e. statres -3) 
DATA_INTEGER(n_groups);      //number groups (i.e., unique combos of LH,stream,downstream,year). Used in psi mlogit backtransform
DATA_INTEGER(n_unique_CH); //number of unique capture occasions
DATA_IMATRIX(CH);           //capture histories (excluding occasion at marking (which we are conditioning on))
DATA_IVECTOR(freq);       //frequency of capture histories
DATA_MATRIX(X_phi);         //fixed effect design matrix for phi
DATA_MATRIX(X_p);         // fixed effect design matrix for psi
DATA_MATRIX(X_psi);         // fixed effect design matrix for psi
DATA_STRUCT(Phi_pim, pim); //index vector of matrices for phi parameter vector for a given Ch x occasion
DATA_STRUCT(p_pim, pim);   //index vector of matrices for p parameter vector for a given Ch x occasion
DATA_IVECTOR(Psi_pim);
DATA_IVECTOR(fix_p_last);
  //~~~~~~~~~~~~~~~~~~~
  // Parameters
  //~~~~~~~~~~~~~~~~~~~
  PARAMETER_VECTOR(beta_phi);  //Phi fixed effect coefficients
  PARAMETER_VECTOR(beta_p);    //p fixed effect coefficients
  PARAMETER_VECTOR(beta_psi);    //psi fixed effect coefficients

  //~~~~~~~~~~~~~~~~~~~
  // Variables
  //~~~~~~~~~~~~~~~~~~~
  
  // Joint negative log-likelihood
  parallel_accumulator<Type> jnll(this);

  // Linear predictors
  vector<Type> eta_phi = X_phi*beta_phi;
  vector<Type> eta_p = X_p*beta_p;
  vector<Type> eta_psi = X_psi*beta_psi;


  // Apply link
  vector<Type> phi=invlogit(eta_phi);
  vector<Type> p=invlogit(eta_p);
  for (int i =0; i<fix_p_last.size(); i++) p(fix_p_last(i))=1; //fix detection at 1 on last occasion
  //** will this crash things if fix_p_last.size() =0?**
  REPORT(phi);
  REPORT(p);
  ////phi inverse multinomial logit
  matrix<Type> psi(n_groups,n_states);
  eta_psi= exp(eta_psi);
  vector<Type> denom = eta_psi.segment(0,n_groups)+eta_psi.segment(n_groups,n_groups)+Type(1);
  psi.col(0)= Type(1)/denom;                      //return after 1 year
  psi.col(1)= eta_psi.segment(0,n_groups)/denom; //return after 2 year
  psi.col(2)= eta_psi.segment(n_groups,n_groups)/denom; //return after 3 year
  REPORT(psi);

  //~~~~~~~~~~~~~~~~~~~
  // Likelihood
  //~~~~~~~~~~~~~~~~~~~

  ////Variables
  vector<Type> pS(4); //state probs: dead, 1, 2, 3
  Type u = 0;         // holds the sum of probs after each occasion
  Type NLL_it=0;      // holds the NLL for each CH
  Type tmp = 0;       // holds the prob of a given state during observation process in upstream migration

  for(int n=0; n<n_unique_CH; n++){ // loop over individual unique capture histories
    pS.setZero(); //initialize at 0,1,0,0 (conditioning at capture)
    pS(1)=Type(1);
    NLL_it=Type(0); //initialize capture history NLL at 0

  //downstream migration
  for(int t=0; t<nDS_OCC; t++){       //loop over downstream occasions (excluding capture occasion)
    //survival process
    pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1)); //prob die or stay dead
    pS(1) *= Type(phi(Phi_pim(0)(n,t))); //prob stay alive
    //observation process
    pS(1) *= Type(p(p_pim(0)(n,t))*CH(n,t)+ (Type(1)-p(p_pim(0)(n,t)))*(Type(1)-CH(n,t))); //prob observation given alive
    pS(0) *= Type(Type(1)-CH(n,t)); //prob observation given dead
    //acculate NLL
    u = pS.sum();  //sum of probs
    pS = pS/u; //normalize probs
    NLL_it  +=log(u);    //accumulate nll
  }

  //ocean occasion
  int t = nDS_OCC;  //set occastion to be ocean occasion
  ////survival process
  pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1)); //prob die or stay dead in ocean
  pS(1) *= Type(phi(Phi_pim(0)(n,t))); //prob survive ocean
  //maturation age process
  pS(2) = pS(1) * psi(Psi_pim(n),1);
  pS(3) = pS(1) * psi(Psi_pim(n),2);
  pS(1) *= psi(Psi_pim(n),0);
  ////observation process
  if(!CH(n,t)){
  pS(1) *= Type(Type(1)-p(p_pim(0)(n,t)));
  pS(2) *= Type(Type(1)-p(p_pim(1)(n,t)));
  pS(3) *= Type(Type(1)-p(p_pim(2)(n,t)));
  } else{
    tmp=pS(CH(n,t))*p(p_pim((CH(n,t)-1))(n,t));
    pS.setZero();
    pS(CH(n,t))=tmp;
  }
  //acculate NLL
  u = pS.sum();  //sum of probs
  pS = pS/u; //normalize probs
  NLL_it  +=log(u);    //accumulate nll
  //end ocean occasion

  //upstream migration
  for(int t=(nDS_OCC+1); t<n_OCC; t++){       //loop over upstream occasions
    ////survival process
    pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1))+
      Type((Type(1)-phi(Phi_pim(1)(n,t)))*pS(2))+
      Type((Type(1)-phi(Phi_pim(2)(n,t)))*pS(3));  // sum(prob vec * 1, 1-phi_1, 1-phi_2, 1-phi_3)
    pS(1) *= Type(phi(Phi_pim(0)(n,t)));                 // sum(prob vec * 0,   phi_1,       0,       0)
    pS(2) *=  Type(phi(Phi_pim(1)(n,t)));                 // sum(prob vec * 0,       0,   phi_2,       0)
    pS(3) *=  Type(phi(Phi_pim(2)(n,t)));                 // sum(prob vec * 0,       0,       0,   phi_3)
    ////observation process
    if(!CH(n,t)){
      pS(1) *=  Type((Type(1)-p(p_pim(0)(n,t))));
      pS(2) *=  Type((Type(1)-p(p_pim(1)(n,t))));
      pS(3) *=  Type((Type(1)-p(p_pim(2)(n,t))));
    }else{
      tmp=pS(CH(n,t))*p(p_pim((CH(n,t)-1))(n,t));
      pS.setZero();
      pS(CH(n,t))=tmp;
    }
    //acculate NLL
  u = pS.sum();  //sum of probs
  pS = pS/u; //normalize probs
  NLL_it  +=log(u);    //accumulate nll
  }
  //multiply the NLL of an individual CH by the frequency of that CH and subtract from total jnll
  jnll-=(NLL_it*freq(n));
  }
  //end of likelihood
  //return jnll
  return(jnll);
}
