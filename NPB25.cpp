#include<RcppArmadillo.h>
#include<R.h>


template<typename eT>
class Datum
{
public:
  
  static const eT pi;       //!< ratio of any circle's circumference to its diameter
  static const eT e;        //!< base of the natural logarithm
  static const eT euler;    //!< Euler's constant, aka Euler-Mascheroni constant
  static const eT gratio;   //!< golden ratio
  static const eT sqrt2;    //!< square root of 2
  static const eT eps;      //!< the difference between 1 and the least value greater than 1 that is representable
  static const eT log_min;  //!< log of the minimum representable value
  static const eT log_max;  //!< log of the maximum representable value
  static const eT nan;      //!< "not a number"
  static const eT inf;      //!< infinity 
  
  // 
  
  static const eT m_u;       //!< atomic mass constant (in kg)
  static const eT N_A;       //!< Avogadro constant
  static const eT k;         //!< Boltzmann constant (in joules per kelvin)
  static const eT k_evk;     //!< Boltzmann constant (in eV/K)
  static const eT a_0;       //!< Bohr radius (in meters)
  static const eT mu_B;      //!< Bohr magneton
  static const eT Z_0;       //!< characteristic impedance of vacuum (in ohms)
  static const eT G_0;       //!< conductance quantum (in siemens)
  static const eT k_e;       //!< Coulomb's constant (in meters per farad)
  static const eT eps_0;     //!< electric constant (in farads per meter)
  static const eT m_e;       //!< electron mass (in kg)
  static const eT eV;        //!< electron volt (in joules)
  static const eT ec;        //!< elementary charge (in coulombs)
  static const eT F;         //!< Faraday constant (in coulombs)
  static const eT alpha;     //!< fine-structure constant
  static const eT alpha_inv; //!< inverse fine-structure constant
  static const eT K_J;       //!< Josephson constant
  static const eT mu_0;      //!< magnetic constant (in henries per meter)
  static const eT phi_0;     //!< magnetic flux quantum (in webers)
  static const eT R;         //!< molar gas constant (in joules per mole kelvin)
  static const eT G;         //!< Newtonian constant of gravitation (in newton square meters per kilogram squared)
  static const eT h;         //!< Planck constant (in joule seconds)
  static const eT h_bar;     //!< Planck constant over 2 pi, aka reduced Planck constant (in joule seconds)
  static const eT m_p;       //!< proton mass (in kg)
  static const eT R_inf;     //!< Rydberg constant (in reciprocal meters)
  static const eT c_0;       //!< speed of light in vacuum (in meters per second)
  static const eT sigma;     //!< Stefan-Boltzmann constant
  static const eT R_k;       //!< von Klitzing constant (in ohms)
  static const eT b;         //!< Wien wavelength displacement law constant
};


// the long lengths of the constants are for future support of "long double"
// and any smart compiler that does high-precision computation at compile-time

template<typename eT> const eT Datum<eT>::pi        = eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);
template<typename eT> const eT Datum<eT>::e         = eT(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274);
template<typename eT> const eT Datum<eT>::euler     = eT(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495);
template<typename eT> const eT Datum<eT>::gratio    = eT(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374);
template<typename eT> const eT Datum<eT>::sqrt2     = eT(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727);



typedef Datum<float>  fdatum;
typedef Datum<double> datum;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec norep( arma::vec A){
  int
  m = A.n_rows,
    n = 0;
  arma::vec
    SA = sort(A),
      B(1),
      aux(1); 
  B(0) = SA(0);
  for(int i=1; i < m ; i++){
    if(B(n) != SA(i)){
      aux(0) = SA(i);
      B.insert_rows(B.n_rows,aux);
      n      = n +1;
    }
  }  
  return B;  
}

// [[Rcpp::export]]
arma::field<arma::vec> dadosrs(arma::mat X ){
  arma::vec
  sind = norep(X.col(0));
  int
    m = sind.n_rows,  /* # de sistemas*/
    k = max(X.col(2)),/* # de tipo de falha (componentes)*/
    a = X.n_rows;     /* #  de linha da base X*/
  arma::field<arma::vec>
    F(m,k);
  arma::vec
    aux(1),
    vux(2);
  aux(0) = 0 ;
  vux(0) = 0 ;
  vux(1) = 0 ;
  for(int j =0; j < m;j++){
    vux(0) = sind(j);  
    for(int q = 1; q<k+1;q++){
      F(j,q-1) = vux;
    }
  }
  for(int l =0 ; l < a;l++){
    for(int j =0; j < m;j++){
      for(int q = 1; q<k+1;q++){
        if((X(l,0) == sind(j)) & ( X(l,2) == q)){
          aux(0) = X(l,1); 
          F(j,q-1)(1) = F(j,q-1)(1)+1;
          F(j,q-1).insert_rows(F(j,q-1).n_rows,aux);
          q = k + 1;
          j = m + 1;
        }
      }
    }  
  }
  return F;
}  

// [[Rcpp::export]]
arma::vec rnp(int n ){
  arma::vec Y(n);
  for(int i=0; i< n;i++)
    Y(i) = R::rnorm(0,1);
  return Y;
}



/*********************************************************************
 *                MCMC - GIBBS                                       * 
 *********************************************************************/
  
  /* GRUPOS*/
  // [[Rcpp::export]]
arma::field<arma::vec> classres(arma::vec res, arma::vec Z ){
  int 
    n    = res.n_rows,
    zast = max(Z),
    verif=0,
    j;
  arma::field<arma::vec>
    F(zast+2);
  arma::vec
    nclass(zast+1);
  arma::mat  
    aux;
  nclass.zeros();
  for ( int i=0; i < n ; i++){
    verif = 0;
    j     = 0;
    while(verif == 0){
      if( j == Z(i)){
        if(nclass(j) > 0){ 
          aux         = res.row(i);
          F(j+1).insert_rows(F(j+1).n_rows,aux);
          nclass(j)   = nclass(j) + 1;  
          verif       = 1;}
        else{
          aux       = res.row(i);
          F(j+1)    = aux;
          nclass(j) = nclass(j) + 1;  
          verif     = 1;}  
        }
      j = j+1;
    }
  }
  F(0) = nclass;
  return F;
}  

/* PARAMETRO */
double condC( arma::vec nu, double a0, double b0 ){
  int
    n = nu.n_rows;
  double
    slnu = sum(log(nu)),
    a1   = n + a0,
    b1   = 1/(b0 - slnu);
  return R::rgamma(a1,b1);
} 

/* DISTRIBUICAO CONDICIONAL NU DADO Y e C*/ 
 // [[Rcpp::export]]
arma::vec condnu( arma::field<arma::vec> nX, double c){
  int 
    T    ,
    zast = nX(0).n_rows;
  arma::vec
    nu(zast),
    nclass    = nX(0),
    cumnclass = cumsum(nclass);
  T    = sum(nclass);
  zast = nclass.n_rows;
  for ( int j=0; j < zast; j++){
    nu(j)  = R::rbeta( nclass(j) + 1, T - cumnclass(j) + c );  
  }  
  return  nu;
}

// [[Rcpp::export]]
arma::vec wnu( arma::vec nu){
  int 
    m    = nu.n_rows;
  arma::vec
    w(m);
  double
    sw;
  w(0)  = nu(0);
  sw    = w(0);
  for ( int j=1; j < m; j++){
    w(j)   = (1- sw) * nu(j);
    sw     = sw + w(j);
  }
  return  w;
}

// [[Rcpp::export]]
arma::vec wnuj(arma::vec w, double  cond, double c){
  double
    aux,
    sumW;
  arma::vec
    wi = w,
    W(1);
    sumW = sum(wi);
  while (sumW < cond){
    aux    = R::rbeta(1,c);
    W(0)   = (1 - sumW) * aux;
    wi.insert_rows(wi.n_rows, W );
    sumW   = sum(wi);
  }
  return  wi;
}


// [[Rcpp::export]]
arma::vec latenteu(arma::vec w, arma::vec Z,int T){
  arma::vec
    u(T);
    u.randu();
  for ( int i=0; i < T; i++){
    u(i) = u(i) * w(Z(i));
  }
  return u;
}


/*HIPERPARAMETROS DA NORMAL - VARIAVEL Y*/

// [[Rcpp::export]]
arma::mat hiperparnormal(  arma::field<arma::vec> F, arma::vec hiperpar, int jast){
  int
  m = F.n_rows-1,
    o;
  double
    u0 = hiperpar(0),
      c0 = hiperpar(1),
      t0 = hiperpar(2),
      d0 = hiperpar(3),
      p,
      sd,
      my;
  arma::vec
    nY = F(0),
      Y;
  arma::mat
    hiper(4,jast);
  for(int i=0;i < jast; i++){
    if( i < m){
      o = nY(i);
    }else{
      o = 0;
    }
    if(o > 0){
      Y          = F(i+1);
      my         = mean(Y);
      sd         = var(Y);
      p          = c0/(c0+o);
      hiper(0,i) = p*u0 + (1-p)*my;
      hiper(1,i) = c0 + o;
      hiper(2,i) = t0 + o*sd + p*o*pow(u0 - my,2);
      hiper(3,i) = d0 + o;
    }else{
      hiper(0,i) = u0;
      hiper(1,i) = c0;
      hiper(2,i) = t0;
      hiper(3,i) = d0;
    }
  }
  return hiper;
}


/* DISTRIBUI플O NORMALGAMMA*/
// [[Rcpp::export]]
arma::vec normalgamma( arma::vec par){
  double
  u0 = par(0),
    c0 = par(1),
    t0 = par(2),
    d0 = par(3),
    va;
  arma::vec 
    am(2);
  am(1) = R::rgamma(d0/2,2/t0);
  am(0) = R::rnorm(0,1);
  va    = 1/am(1);
  am(0) = u0 + sqrt(va/c0)*am(0);
  return am;
}

arma::mat hipY(arma::mat hiper){
  int
  m = hiper.n_cols;
  arma::mat
    A(2,m);
  for(int i=0; i <m; i++){
    A.col(i) = normalgamma( hiper.col(i));
  }
  return A;
}

/* 4 Sampling Y */
// [[Rcpp::export]]
double ldensnormal(double Z, arma::vec hiperpar ){
  double
  mu  = hiperpar(0),
    var = 1/hiperpar(1);
  return  -0.5*log(var)  - 0.5*pow(Z - mu,2)/var;
}

// [[Rcpp::export]]
arma::vec condY( arma::vec X, arma::vec w, arma::vec u, arma::mat nusig ){
  int
     j   = 0,
    sel = 0,
    T   = X.n_rows;
  double
    sw = sum(w),
      fw,  
      psim,
      verif;
  arma::vec
    paux(1),
    p,
    ipaux(1),
    ip,
    Z(T);
  for ( int i=0; i < T; i++){
    j     = 0;
    p     = paux;
    ip    = ipaux; 
    verif = 0;
    fw    = sw;
    
    while(fw    > u(i)){
      if(w(j) > u(i) ){ 
        if( verif != 0 ){
          paux(0) = ldensnormal( X(i), nusig.col(j) ) ;
          p.insert_rows(p.n_rows,paux);
          ipaux(0) = j;
          ip.insert_rows(ip.n_rows,ipaux);
        }else{
          p(0)  = ldensnormal( X(i), nusig.col(j) ) ;
          ip(0) = j;
          verif = 1;
        }
      }
      fw = fw -w(j);
      j  = j + 1;
    }
    sel  = 0;
    j    = 0;
    p    = p - max(p);
    p    = exp(p);
    p    = p/sum(p);
    p    = cumsum(p);
    psim = R::runif(0, 1 );
    while ( sel == 0){
      if( psim < p(j)){Z(i) = ip(j); sel = 1;}
      else{j = j + 1; }   
    }
  }
  return Z;
}  



// [[Rcpp::export]]
arma::mat amlamb( arma::vec Y, arma::mat nusig, int m ){
  arma::mat
  lambda(2,m);
  for ( int i=0; i < m; i++){
    lambda.col(i) = nusig.col(Y(i));
  }
  return lambda;
}    


// [[Rcpp::export]]
arma::vec simY(int n,double C ){
  int
  aux,
  controle ;
  double
    nu,  
    prodnu = 1,
    umax;
  arma::vec
    vaux(1),
    w(1),
    p,
    u(n),
    Z(n);
  u.randu();
  umax   = max(u);
  nu     = R::rbeta(1,C);
  prodnu = prodnu * (1 - nu); 
  w(0)   = nu;
  while( sum(w) < umax ){
    nu      = R::rbeta(1,C);
    vaux(0) = prodnu * nu;
    w.insert_rows(w.n_rows,vaux);
    prodnu  = prodnu * (1 - nu);
  }
  for(int i=0; i < n; i++){
    p        = cumsum(w);
    controle = 0;
    aux      = 0;
    while(controle ==0){
      if( u(i) < p(aux)){Z(i) = aux;controle = 1;}
      else{aux = aux + 1; }   
    }
  }   
  return Z;      
}  


// [[Rcpp::export]]
arma::vec vecgamma(int n, double a, double  b){
  double
    ib = 1/b;
  arma::vec
    x(n);
  for(int i = 0; i < n; i++){
    x(i) = R::rgamma(a,ib);
  }
  return x;
}


arma::vec vgamma(int k, double a , double b ){
  arma::vec
    A(k);
  for(int i=0 ; i <k;i++){
    A(i) = R::rgamma(a,1/b);
  }
  return A;
}

/*
 * ESTIMADOR DE MAXIMA VEROSSIMILHANCA DO MODELO SIMPLES
 */

// [[Rcpp::export]]
arma::field<arma::mat> MVSM(arma::mat X, double T, int r){
  double
  logT = log(T);
  arma::field<arma::vec> 
    F = dadosrs(X);
  int
    m    = F.n_rows,
    m1   = m + r,  
    k    = F.n_cols;
  arma::vec
    nc = arma::zeros(k),   /* # DE DEFEITOS DE CADA COMPONENTE*/ 
    nj = arma::zeros(m),   /* # DE DEFEITOS DE CADA SISTE*/
    lt = arma::zeros(k),   /* SOMA DO LOG DO INSTANTES DE DEFEITO DE CADA COMPONENTE*/
    alpha(k),
    beta(k),
    lijt;
  arma::field<arma::mat>
    est(5);
  arma::mat
    njq = arma::zeros(m,k);
  for(int i =0; i < m;i++){
    for(int j =0; j < k;j++){
      njq(i,j) = F(i,j)(1);
      if( F(i,j)(1) > 0 ){
        lijt   = F(i,j);
        lijt   = lijt.rows(2,lijt.n_rows -1);  
        nc(j)  = nc(j) + F(i,j)(1);
        lt(j)  = lt(j) + sum(logT - log(lijt) );
        nj(i)  = nj(i) + F(i,j)(1); 
      }  
    }
  }  
  for(int j =0; j < k;j++){
    beta(j)   = nc(j)/lt(j);
    alpha(j)  = logT - (log(nc(j)) - log(m1))/beta(j); /*logaritmo*/
    alpha(j)  = exp(alpha(j));
  }
  est(0) = alpha;
  est(1) = beta;
  est(2) = nc;
  est(3) = nj;
  est(4) = njq;
  return est;
}



// [[Rcpp::export]]
double lvfragm( arma::vec Z, arma::vec n, int k){
  arma::vec
    W      = log(Z),
    lv     =  n % W.rows(0,n.n_rows-1) ;
  return  sum(lv);
}

// [[Rcpp::export]]
arma::vec fsfragm( arma::vec Z,  arma::vec n, int k){
  arma::vec
    aux,
    fs  = arma::zeros(Z.n_rows),
    IZ  = 1/Z.rows(0,n.n_rows-1);
  fs.rows(0,n.n_rows-1) = n % IZ;
  return fs ;
}  


// [[Rcpp::export]]
double prZ(arma::vec Z, arma::mat hiperpar){
  arma::vec
    W    = log(Z),
    mu   = hiperpar.col(0),
    pre  = hiperpar.col(1),
    Y    = W - mu;
  arma::mat
    delta = -0.5* Y.t() * (pre % Y);
  return  - sum(W) + 0.5*sum(log(pre)) + as_scalar(delta);
}

// [[Rcpp::export]]
arma::vec dprZ(arma::vec Z, arma::mat hiperpar){
  arma::vec
    W    = log(Z),
    IZ   = 1/Z,
    mu   = hiperpar.col(0),
    pre  = hiperpar.col(1),
    Y    = W - mu,
    ddelta = - pre % Y;
  return  -IZ +  ddelta % IZ ;
}

// [[Rcpp::export]]
arma::vec reparz(arma::vec y){
  int
    m = y.n_rows;
  arma::vec
    a, 
    z,
    aux(m);
  for(int i=0;i < m;i++ ){
    aux(i) = - log(m - i);
  }
  a  =  y + aux;
  z  = 1/(1+exp(-a));
  return z;
}  
  
// [[Rcpp::export]]
arma::vec dreparz(arma::vec y){
  arma::vec
    z  = reparz( y);
  return z % (1- z);
}


// [[Rcpp::export]]
arma::vec reparx(arma::vec z){
  int
    m = z.n_rows;
  double
    sx = 0;
  arma::vec
    x(m);  
  for(int i=0;i < m;i++ ){
    x(i) = (1 - sx) *z(i);
    sx   = sx + x(i);
  }
  return x;
}


// [[Rcpp::export]]
arma::vec dreparx( arma::vec x){
  int
    m = x.n_rows;
  double
    sx = 0;
  arma::vec
    dz(m);  
  for(int i=0;i < m;i++ ){
    dz(i) =  (1 - sx);
    sx   = sx + x(i);
  }
  return dz;
}



// [[Rcpp::export]]
arma::vec drepardx(arma::vec dx , arma::vec x, arma::vec z){
  int
    m = x.n_rows;
  double
    aux,
    dev;
  arma::vec
    dxz = dreparx( x),
    dxv,
    dxx(m);  
  for(int i=0;i < m - 1;i++ ){
    dxv    = arma::zeros(m);
    dxv(i) = dxz(i);
    aux    = dxv(i);
    dev    = dx(i) * dxv(i); 
    for(int j=i+1;j < m ;j++ ){
      dxv(j) = - z(j) *aux ; 
      aux    = aux + dxv(j);
      dev    = dev + dx(j) * dxv(j);  
    }  
    dxx(i) =  dev;
  }
  dxx(m-1) = dx(m-1) * dxz(m-1) ;
  return dxx;
}


// [[Rcpp::export]]
double jacobum(arma::vec x, arma::vec z, int m){
  arma::vec
    dz = z % (1-z);
  double 
    sw = 0,
    lp = sum(log(dz));
  for(int i=0;i < m ;i++ ){
    sw = sw + x(i);
    lp = lp + log(1 - sw);
  }
  return lp;
}
  

  // [[Rcpp::export]]
arma::vec djacobum(arma::vec x, arma::vec z, int m){
  arma::vec
    dz = z % (1-z),
    fs,
    fx = arma::zeros(m),
    dx = dreparx(x),
    dxx;
  double 
    sw = 0;
  for(int i=0;i < m ;i++ ){
    sw           = sw + x(i);
    for(int j = 0; j < i + 1; j++ ){
      fx(j) =  fx(j)  - 1/(1-sw)  ;
    }  
  }
  dxx           = drepardx(fx, x, z);
  fs            = (1/z - 1/(1-z) + dxx) % dz;
  return fs;
  }
  
  
  
// [[Rcpp::export]]
arma::vec fsfragmY( arma::vec y,  arma::vec alpha, arma::vec n, int k, arma::mat hiperpar){
  int
    m = y.n_rows + 1;
  arma::vec
    Um = arma::ones(m),
    ey = exp(-y),
    z  = reparz(y) ,
    x  = reparx(z),
    dz = z % (1-z),
    dx = dreparx(x),
    Z(m),
    fx,
    fZ,
    f1,
    dxx,
    fs,
    djacob        = djacobum(x,z ,m-1);
    Z.rows(0,m-2) = m*x;
    Z(m-1)        = m  - sum(Z.rows(0,m-2));
    fZ            = dprZ( Z, hiperpar) + fsfragm( Z, n, k);
    f1            = fZ.rows(0,m-2)- fZ(m-1);  
    fx            = m *f1;
    dxx           = drepardx(fx, x, z);
    fs            = dxx % dz + djacob;
  return  fs;
}


/*############################################################################
#  LEAPFROG                                                                  #
############################################################################*/

// [[Rcpp::export]]
arma::mat LeapfrogfragmY (arma::mat thetar, arma::vec alpha,  double epsilon, arma::vec n, arma::mat hiperpar, int k ){
  arma::mat 
  thetan = thetar;
  thetan.col(1)  = thetar.col(1) + epsilon  * fsfragmY( thetar.col(0),  alpha, n, k, hiperpar);  
  thetan.col(0)  = thetar.col(0) + epsilon  * thetan.col(1);
  return thetan;
}


/*############################################################################
#            A FUN플O DENSIDADE DE PROBABILIDADE                             #  
############################################################################*/

// [[Rcpp::export]]
double fdensfragmY(arma::mat thetar, arma::vec alpha, arma::vec n, arma::mat hiperpar,int k){
  arma::vec
    y   = thetar.col(0), 
    r   = thetar.col(1);
  int
    m = y.n_rows + 1;
  arma::vec
    ey  = exp(-y),
    z  = reparz(y) ,
    x  = reparx(z),
    Z(m);
  Z.rows(0,m-2) = m*x;
  Z(m-1)        = m  - sum(Z.rows(0,m-2));
  arma::mat
    delta =  - 0.5* r.t()  * r;
  double
    jacob = jacobum(x, z, m-1),
    lp    = lvfragm(  Z, n, k)     + prZ(Z, hiperpar) 
          + jacob 
          + as_scalar(delta);
    return lp;     
}


/*#####################################################
#            HMC                                     #
#####################################################*/

// [[Rcpp::export]]
arma::vec amHMCfragmY( arma::vec theta_current, arma::vec alpha, arma::vec n, double epsilon, int LF, arma::mat hiperpar, int k, int m){
  double 
    ep = 0.5 * epsilon,
    p,
    H_current,
    H_prop,
    limit = - exp(999);
  arma::mat 
    thetan(m-1 ,2);
  arma::vec
    ace(1), 
    theta(m-1); 
  ace(0) = 0;
  theta         = theta_current;
  thetan.col(0) = theta;
  thetan.col(1) = rnp(m-1);
  H_current     = fdensfragmY(thetan, alpha, n, hiperpar, k);
  thetan.col(1) = thetan.col(1) + ep* fsfragmY( thetan.col(0),  alpha, n, k, hiperpar);
  thetan.col(0) = thetan.col(0) + epsilon  * thetan.col(1);
  for (int j=0; j < LF - 1; j++){
    thetan  = LeapfrogfragmY (thetan,  alpha, epsilon, n, hiperpar, k );
  }
  thetan.col(1) = thetan.col(1) + ep* fsfragmY( thetan.col(0),  alpha, n, k, hiperpar);
  H_prop        = fdensfragmY(thetan, alpha, n, hiperpar, k);
  p             = R::runif(0, 1 );
  if(thetan.has_nan() == 0){
    if ( (log(p) < H_prop - H_current) & ( H_prop > limit) ){
      theta       = thetan.col(0);
      ace(0)      = 1;
    }
  }  
  theta.insert_rows(theta.n_rows,ace);
  return theta ;   
}


// [[Rcpp::export]]
arma::vec Ynodup(arma::vec A){
  int
  m  = A.n_rows;
  arma::vec
    SA = sort(A),
      B(1),
      aux(1);
  aux(0) = SA(0);
  B(0)   = aux(0); 
  for(int i=1; i < m;i++){
    if( aux(0) != SA(i)){
      aux(0) = SA(i);
      B.insert_rows(B.n_rows,aux);
    } 
  }
  return B;
} 

// [[Rcpp::export]]
arma::vec Yorder(arma::vec A){
  int
  m  = A.n_rows,
    n,
    ref;
  arma::vec
    SY = Ynodup(A),
      B(m); 
  for(int i=0; i < m;i++){
    ref = 0;
    n  = 0 ;
    while(ref == 0){
      if(A(i) == SY(n)){
        B(i) = n;
        ref  = 1;
      }else{n = n + 1;} 
    }
  }
  return B;
} 

/* MCMC */
// [[Rcpp::export]]
arma::field<arma::mat> GHMCRSMNB( arma::mat X, double T, int r, int warmup, int iter, double epsilon , int LF, arma::vec parhiperc, arma::vec parnormal){ 
  arma::field<arma::vec> 
    F = dadosrs(X);
  int
    mc   = F.n_rows,  
    m    = mc  + r,
    k    = F.n_cols,
    jast,
    ref  = 0,
    SS   = warmup + iter;
  double
    im    = m,
    cont  = 0.1*SS,
    pcont = cont;
  arma::vec
    para(k),
    pary(m-1),
    parZ(m),
    nu,
    w,
    u,
    Y,
    z,
    x,
    aux,
    SUZ(iter),
    LZ(m),
    pac =arma::zeros(1),
    AC(iter),
    nj,
    nq;
  double
    C = R::rgamma(1,1),
    cond;
  arma::mat
    njq = arma::zeros(m,k),
    alpha(k,iter),
    beta(k,iter),
    W(m,iter),
    hipnu,
    nuY,
    nusig,
    lv(m,iter),
    auxZ;
  arma::field<arma::vec> 
    nY;
  arma::field<arma::mat>
    mcmc(iter+2,3),
    estclass =  MVSM(X, T, r);
  pary.zeros();
  nq               = estclass(2);
  if(r > 0) {
    nj     = estclass(3);
    
  }else{
    nj    = estclass(3);
    nj    = nj.rows(0,mc-2);
  }
  im               = 1/im;
  z                = reparz(pary);
  x                = reparx(z);
  parZ.rows(0,m-2) = m*x;
  parZ(m-1)        = m  - sum(parZ.rows(0,m-2));
  LZ               = log(parZ);
  Y                = simY( m, C );
  Rcpp::Rcout  << "............................... " << ref << "%" << "(warmup)" <<  std::endl;
  ref  = ref + 10;
  for(int t = 0; t < warmup;t++){
    Y            = Yorder(Y); /* CLASSIFICA플O DE CADA CLUSTER */
    nY           = classres(LZ, Y);
    /* SAMPLING - C*/
    C            = condC(nu,  parhiperc(0), 1/parhiperc(1));
    /* SAMPLING - nu */
    nu           = condnu(nY, C);
    w            = wnu( nu);
    /* SAMPLING - u */
    u            = latenteu(w, Y, m);
    cond         = 1 - min(u);
    /* SAMPLING - w */
    w            = wnuj( w, cond, C);
    jast         = w.n_rows;
    hipnu        = hiperparnormal( nY,  parnormal,jast);
    nusig        = hipY(hipnu);
    /* SAMPLING - Y */
    Y        = condY( LZ, w, u, nusig );
    /* SAMPLING - Z*/
    nuY              = amlamb( Y, nusig, m );
    pary             = amHMCfragmY( pary, para,  nj,  epsilon, LF, nuY.t(), k, m);
    pary             = pary.rows(0,m-2);
    z                = reparz(pary) ;
    x                = reparx(z);
    parZ.rows(0,m-2) = m*x;
    parZ(m-1)        = m  - sum(parZ.rows(0,m-2));
    LZ               = log(parZ);
    if( t > cont){
      Rcpp::Rcout  << "..............................." << ref << "%" <<  "(warmup)" <<std::endl;
      cont = cont + pcont;
      ref  = ref + 10;
    }
  }
  for(int t = 0; t < iter;t++){
    Y            = Yorder(Y); /* CLASSIFICA플O DE CADA CLUSTER */
    nY           = classres(LZ, Y);
    /* SAMPLING - C*/
    C            = condC(nu,  parhiperc(0), 1/parhiperc(1));
    AC(t)        = C;
    /* SAMPLING - nu */
    nu           = condnu(nY, C);
    w            = wnu( nu);
    /* SAMPLING - u */
    u            = latenteu(w, Y, m);
    cond         = 1 - min(u);
    /* SAMPLING - w */
    w            = wnuj( w, cond, C);
    jast         = w.n_rows;
    hipnu        = hiperparnormal( nY,  parnormal,jast);
    nusig        = hipY(hipnu);
    mcmc(t+2,0)  = w;
    mcmc(t+2,1)  = nusig;
    /* SAMPLING - Y */
    Y          = condY( LZ, w, u, nusig );
    /* SAMPLING - Z*/
    for(int j = 0; j < k; j++){
      para(j)  = R::rgamma(nq(j), im); 
    }
    nuY              = amlamb( Y, nusig, m);
    pary             = amHMCfragmY( pary, para,  nj,  epsilon, LF, nuY.t(), k, m);
    pac(0)           = pac(0) + pary(m-1);  
    pary             = pary.rows(0,m-2);
    z                = reparz(pary) ;
    x                = reparx(z);
    parZ.rows(0,m-2) = m*x;
    parZ(m-1)        = m  - sum(parZ.rows(0,m-2));
    LZ               = log(parZ);
    W.col(t)         = parZ;
    if( t + warmup > cont){
      Rcpp::Rcout  << "..............................." << ref << "%" <<  "( iter )" << std::endl;
      cont = cont + pcont;
      ref  = ref + 10;
    }
  }
  mcmc(0,0)    = W   ;
  mcmc(0,1)    = pac/iter;
  mcmc(0,2)    = AC      ;
  mcmc(1,1)    = Y;
  return mcmc;
}


// [[Rcpp::export]]
arma::vec vlnormal(arma::vec x, double mu, double sig){
  int
  a = x.n_rows;
  arma::vec 
    fdens =arma::zeros(a);
  for(int i=0; i < a; i++){
    fdens(i) = R::dlnorm(x(i), mu,sig,0);
  }
  return fdens;
}

// [[Rcpp::export]]
arma::mat fdensfraZ(arma::vec x, arma::field<arma::mat> W, arma::field<arma::mat> parN,  int R, arma::vec parnormal){ 
  int
  m = W.n_rows - 2,
    a = x.n_rows,
    find,
    nw,
    l;
  double
    ru,
    mu,
    sig;
  arma::vec
    w,
    u(R),
    f1,
    f,
    aux;
  arma::mat
    nusig,
    fm(a,m);
  for(int i=0; i < m; i++){
    u.randu();
    u     = sort(u);
    w     = W(i+2);
    w     = cumsum(w);
    nw    = w.n_rows;
    nusig = parN(i+2);
    l     = 0;
    f   = arma::zeros(a);
    for(int j=0; j < R; j++){
      ru   = u(j);
      find = 0 ;
      while( find == 0){
        if( l >= nw){
          aux  = normalgamma( parnormal );
          mu   = aux(0);
          sig  = 1/aux(1);
          f1   = vlnormal( x,  mu, sig);
          f    = f + f1/R;
          find = 1;
        }else{if( ru < w(l)){
          mu   = nusig(0,l);    
          sig  = 1/nusig(1,l);
          f1   = vlnormal( x,  mu, sig);
          f    = f + f1/R;
          find = 1;
        }else{l =  l + 1;}}}
    }
    fm.col(i) = f;
  }
  return  mean(fm,1);
}   

// [[Rcpp::export]]    
arma::vec mcent(arma::vec mon){
  arma::vec
  mc(4);
  mc(0) = mon(0),
    mc(1) = mon(1) - pow(mon(0),2),
    mc(2) = mon(2) - 3*mon(0)*mon(1) + 2*pow(mon(0),3),
    mc(3) = mon(3) - 4*mon(0)*mon(2) + 6*pow(mon(0),2)*mon(1) - 3*pow(mon(0),4);  
  return mc;
}   


// [[Rcpp::export]]   
arma::field<arma::mat> estBZ( arma::field<arma::mat> W, arma::field<arma::mat> parN,  int R, arma::vec parnormal){ 
  int
    m = W.n_rows - 2,
    find,
    nw,
    l;
  double
    ru,
    mu,
    sig;
  arma::vec
    w,
    u(R),
    mw0(4),
    mz0(4),
    mw1(4),
    mz1(4),
    aux;
  arma::mat
    nusig,
    mw(4,m),
    mz(4,m);
  arma::field<arma::mat>
    estBZ(8);
  for(int i=0; i < m; i++){
    u.randu();
    u     = sort(u);
    w     = W(i+2);
    w     = cumsum(w);
    nw    = w.n_rows;
    nusig = parN(i+2);
    l     = 0;
    mw0.zeros(); 
    mz0.zeros();
    mw1.zeros(); 
    mz1.zeros();
    for(int j=0; j < R; j++){
      ru   = u(j);
      find = 0 ;
      while( find == 0){
        if( l >= nw){
          aux    = normalgamma( parnormal );
          mu     = aux(0);
          sig    = 1/aux(1);
          mw0(0) = mu;
          mw0(1) = pow(mu,2) + sig;
          mw0(2) = pow(mu,3) + 3*mu*sig;
          mw0(3) = pow(mu,4) + 6*pow(mu,2)*sig + 3*sig;
          mz0(0) = exp(   mu + 0.5*sig);
          mz0(1) = exp( 2*mu + 2  *sig);
          mz0(2) = exp( 3*mu + 4.5*sig);
          mz0(3) = exp( 4*mu + 8  *sig);
          mw1    = mw1 + mw0/R; 
          mz1    = mz1 + mz0/R; 
          find = 1;
        }else{if( ru < w(l)){
          mu     = nusig(0,l);    
          sig    = 1/nusig(1,l);
          mw0(0) = mu;
          mw0(1) = pow(mu,2) + sig;
          mw0(2) = pow(mu,3) + 3*mu*sig;
          mw0(3) = pow(mu,4) + 6*pow(mu,2)*sig + 3*sig;
          mz0(0) = exp(   mu + 0.5*sig);
          mz0(1) = exp( 2*mu + 2  *sig);
          mz0(2) = exp( 3*mu + 4.5*sig);
          mz0(3) = exp( 4*mu + 8  *sig);
          mw1    = mw1 + mw0/R; 
          mz1    = mz1 + mz0/R;
          find = 1;
        }else{l =  l + 1;}}}
    }
    mw.col(i) = mw1;
    mz.col(i) = mz1;
  }
  estBZ(0)    = mean(mw,1);
  estBZ(1)    = mean(mz,1);
  estBZ(2)    = mw;
  estBZ(3)    = mz;
  estBZ(4)    = w;
  estBZ(5)    = u;
  estBZ(6)    = mz0;
  estBZ(7)    = arma::zeros(2);
  return  estBZ;  
}   

