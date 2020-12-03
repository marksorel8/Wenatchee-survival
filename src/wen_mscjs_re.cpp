#include <TMB.hpp>

 
//Multistate cormacck jolly seber model for downstream migrating salmon in the Columbia River
// 
// Copyright (C) 2020  Mark Sorel
// 
// This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU Affero General Public License as
//   published by the Free Software Foundation, either version 3 of the
//   License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Affero General Public License for more details.
//   
//   You should have received a copy of the GNU Affero General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//     
//     I can be contacted at marks6@uw.edu or at:
//       Mark Sorel
//       1122 NA Boat Street,
//       Seattle, WA 98105
// 


// This program was inspired by packages marked by Jeff Laake, Devin Johnson, and Paul Conn 
// <https://cran.r-project.org/web/packages/marked/index.html>  and
// package glmmTMB by Mollie Brooks et al. <https://cran.r-project.org/web/packages/glmmTMB/index.html>
// Code chunks were copied directly from the glmmTMB source code with permission form
// some glmmTMB coauthors. 



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Functions copied from glmmTMB for calculating random effect probabilities
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//valid covariance structures
enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3,
  ou_covstruct   = 4,
  exp_covstruct = 5,
  gau_covstruct = 6,
  mat_covstruct = 7,
  toep_covstruct = 8
};

//defines elements of list data structure
template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};


//translates r list data structure to C/TMB list data structure. 
//Returns a vector of list, where each  list is for one random effect "component" i.e. inside one parenthesid e.g.(LH|Year)
template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!isNull(t)){
        RObjectTestExpectedType(t, &isNumeric, "times");
        (*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!isNull(d)){
        RObjectTestExpectedType(d, &isMatrix, "dist");
        (*this)(i).dist = asMatrix<Type>(d);
      }
    }
  }
};


//function that calculates the probability of random effects for many different random effects structures.
//Returns negative log prob of random effects for a given random effect compnenet e.g. (LH|year)
template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      if (do_simulate) {
        U.col(i) = rnorm(Type(0), sd);
      }
    }
    term.sd = sd; // For report
  }
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == cs_covstruct){
    // case: cs_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    Type corr_transf = theta(n);
    vector<Type> sd = exp(logsd);
    Type a = Type(1) / (Type(n) - Type(1));
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) : rho);
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == toep_covstruct){
    // case: toep_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> sd = exp(logsd);
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) :
                       parms( (i > j ? i-j : j-i) - 1 ) );
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == ar1_covstruct){
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    int n = term.blockSize;
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
        ans -= dnorm(U(i, j), phi * U(i-1, j), sd * sqrt(1 - phi*phi), true);
        if (do_simulate) {
          U(i, j) = rnorm( phi * U(i-1, j), sd * sqrt(1 - phi*phi) );
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
        term.sd(i) = sd;
        for(int j=0; j<n; j++){
          term.corr(i,j) = pow(phi, abs(i-j));
        }
      }
    }
  }
  else if (term.blockCode == ou_covstruct){
    // case: ou_covstruct
    //  * NOTE: this is the continuous time version of ar1.
    //          One-step correlation must be non-negative
    //  * NOTE: 'times' assumed sorted !
    int n = term.times.size();
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
        Type rho = exp(-exp(corr_transf) * (term.times(i) - term.times(i-1)));
        ans -= dnorm(U(i, j), rho * U(i-1, j), sd * sqrt(1 - rho*rho), true);
        if (do_simulate) {
          U(i, j) = rnorm( rho * U(i-1, j), sd * sqrt(1 - rho*rho));
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
        term.sd(i) = sd;
        for(int j=0; j<n; j++){
          term.corr(i,j) =
            exp(-exp(corr_transf) * CppAD::abs(term.times(i) - term.times(j)));
        }
      }
    }
  }
  // Spatial correlation structures
  else if (term.blockCode == exp_covstruct ||
           term.blockCode == gau_covstruct ||
           term.blockCode == mat_covstruct){
    int n = term.blockSize;
    matrix<Type> dist = term.dist;
    if(! ( dist.cols() == n && dist.rows() == n ) )
      error ("Dimension of distance matrix must equal blocksize.");
    // First parameter is sd
    Type sd = exp( theta(0) );
    // Setup correlation matrix
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        switch (term.blockCode) {
        case exp_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -dist(i,j) * exp(-theta(1)) ) );
          break;
        case gau_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -pow(dist(i,j),2) * exp(-2. * theta(1)) ) );
          break;
        case mat_covstruct:
          corr(i,j) = (i==j ? Type(1) : matern( dist(i,j),
                       exp(theta(1)) /* range */,
                       exp(theta(2)) /* smoothness */) );
          break;
        default:
          error("Not implemented");
        }
      }
    }
    density::MVNORM_t<Type> nldens(corr);
    density::SCALE_t<density::MVNORM_t<Type> > scnldens = density::SCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = corr;   // For report
    term.sd.resize(n);  // For report
    term.sd.fill(sd);
  }
  else error("covStruct not implemented!");
  return ans;
}


//function that creats the structures and call termwise_nll for all random effects. 
//Returns negative log prob of random effects.
template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
                  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}


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
DATA_INTEGER(Est_J_time); // capture occasion of estuary Juvenile
//CH data
DATA_IMATRIX(CH);           //capture histories (excluding occasion at marking (which we are conditioning on))
DATA_IVECTOR(freq);       //frequency of capture histories
//design matrices fixed effects
DATA_MATRIX(X_phi);         //fixed effect design matrix for phi
DATA_MATRIX(X_p);         // fixed effect design matrix for psi
DATA_MATRIX(X_psi);         // fixed effect design matrix for psi
//design matrices random effects
DATA_SPARSE_MATRIX(Z_phi); //random effect design matrix for phi
//PIMS
DATA_STRUCT(Phi_pim, pim); //index vector of matrices for phi parameter vector for a given Ch x occasion
DATA_STRUCT(p_pim, pim);   //index vector of matrices for p parameter vector for a given Ch x occasion
DATA_IVECTOR(Psi_pim);
//DATA_IVECTOR(fix_p_last);  //indices of last occasion for p, to fix at 1
// Covariance structures 
DATA_STRUCT(phi_terms, terms_t);//  Covariance structure for the Phi model



  //~~~~~~~~~~~~~~~~~~~
  // Parameters
  //~~~~~~~~~~~~~~~~~~~
  //fixed effects
  PARAMETER_VECTOR(beta_phi);  //Phi fixed effect coefficients
  PARAMETER_VECTOR(beta_p);    //p fixed effect coefficients
  PARAMETER_VECTOR(beta_psi);    //psi fixed effect coefficients
  //random effects
  PARAMETER_VECTOR(b_phi);      //phi random effects
  // Joint vector of covariance parameters
  PARAMETER_VECTOR(theta_phi);   //phi
  
  //~~~~~~~~~~~~~~~~~~~
  // Variables
  //~~~~~~~~~~~~~~~~~~~
  
  // Joint negative log-likelihood
  parallel_accumulator<Type> jnll(this);

  // Linear predictors
  //// Fixed component
  vector<Type> eta_phi = X_phi*beta_phi;
  vector<Type> eta_p = X_p*beta_p;
  vector<Type> eta_psi = X_psi*beta_psi;
  //// Random component
   eta_phi += Z_phi*b_phi;
  // eta_p = X_p*beta_p;
  // eta_psi = X_psi*beta_psi;

  // Apply link
  vector<Type> phi=invlogit(eta_phi);
  vector<Type> p=invlogit(eta_p);
  // for (int i =0; i<fix_p_last.size(); i++) p(fix_p_last(i))=1; //fix detection at 1 on last occasion
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
  
  // Random effects
  jnll += allterms_nll(b_phi, theta_phi, phi_terms, this->do_simulate);
  
  

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
    if(t!=Est_J_time){                // fix survival at 1 between Bon_J and Est_J
    pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1)); //prob die or stay dead
    pS(1) *= Type(phi(Phi_pim(0)(n,t))); //prob stay alive
    }
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
  
  
  for(int t=(nDS_OCC+1); t<n_OCC; t++){       //loop over upstream occasions  
  
  ////observation process at t-1 (Obs_t below), because I'm going to fix the detection prob at 1 for the last occasion after this loop
  int Obs_t=t-1;
  if(!CH(n,t-1)){
  pS(1) *= Type(Type(1)-p(p_pim(0)(n,Obs_t)));
  pS(2) *= Type(Type(1)-p(p_pim(1)(n,Obs_t)));
  pS(3) *= Type(Type(1)-p(p_pim(2)(n,Obs_t)));
  } else{
    tmp=pS(CH(n,Obs_t))*p(p_pim((CH(n,Obs_t)-1))(n,Obs_t));
    pS.setZero();
    pS(CH(n,Obs_t))=tmp;
  }
  //accumlate NLL
  u = pS.sum();  //sum of probs
  pS = pS/u; //normalize probs
  NLL_it  +=log(u);    //accumulate nll
  //end ocean occasion

  //upstream migration

    ////survival process at time t
    pS(0) += Type((Type(1)-phi(Phi_pim(0)(n,t)))*pS(1))+
      Type((Type(1)-phi(Phi_pim(1)(n,t)))*pS(2))+
      Type((Type(1)-phi(Phi_pim(2)(n,t)))*pS(3));  // sum(prob vec * 1, 1-phi_1, 1-phi_2, 1-phi_3)
    pS(1) *= Type(phi(Phi_pim(0)(n,t)));                 // sum(prob vec * 0,   phi_1,       0,       0)
    pS(2) *=  Type(phi(Phi_pim(1)(n,t)));                 // sum(prob vec * 0,       0,   phi_2,       0)
    pS(3) *=  Type(phi(Phi_pim(2)(n,t)));                 // sum(prob vec * 0,       0,       0,   phi_3)
  }
    
    ////observation process at final time assuming detection probability is 1
    if(!CH(n,(n_OCC-1))){
      pS(1) =  Type(0);
      pS(2) =  Type(0);
      pS(3) =  Type(0);
    }else{
      tmp=pS(CH(n,(n_OCC-1)));
      pS.setZero();
      pS(CH(n,(n_OCC-1)))=tmp;
    }
    //accmulate NLL
  u = pS.sum();  //sum of probs
  pS = pS/u; //normalize probs
  NLL_it  +=log(u);    //accumulate nll
  //end observation process at final time
  
  //multiply the NLL of an individual CH by the frequency of that CH and subtract from total jnll
  jnll-=(NLL_it*freq(n));
  }
  //end of likelihood
  //return jnll
  return(jnll);
}
