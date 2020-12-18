#include <TMB.hpp>

//cormacck jolly seber model for downstream migrating salmon in the Columbia River
// Wrriten my Mark Sorel, but inspired by packages marked (J. Laake et.al) and
//glmmTMB (M. Brooks et al.).
//considerable code taken from glmmTMB currently without permission permission

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Functions copied from glmmTMB for calculating random effect probabilities
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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