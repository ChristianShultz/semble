// main.cc -
//
// Wednesday, April 11 2012
//

#include "semble/semble_semble.h"
#include "itpp/itbase.h"
#include <vector>
#include <iostream>
#include <string>

using namespace SEMBLE;

template<typename T> 
bool test_equiv(const itpp::Vec<T> &v, 
    const SEMBLE::SembleVector<T> &S,
    const double thresh=1e-6)
{
  itpp::Vec<T> zero(v); 
  zero.zeros(); 

  for ( int i = 0 ; i < S.bins(); ++i)
    if( itpp::round_to_zero(v - S[i],thresh) != zero)
    {
      std::cout << "****\n\n"
        << __func__ << " bin = " << i 
        << std::endl;
      std::cout << v << std::endl; 
      std::cout << "\n" << S[i] << std::endl;
      std::cout << "\n" 
        << itpp::round_to_zero(v - S[i], thresh) << std::endl; 
      return false; 
    }
  return true; 
}


template<typename T> 
bool test_equiv(const itpp::Mat<T> &v, 
    const SEMBLE::SembleMatrix<T> &S,
    const double thresh=1e-6)
{
  itpp::Mat<T> zero(v); 
  zero.zeros(); 

  for ( int i = 0 ; i < S.bins(); ++i)
    if( itpp::round_to_zero(v - S[i],thresh) != zero)
    {
      std::cout << "****\n\n"
        << __func__ << " bin = " << i 
        << std::endl;
      std::cout << v << std::endl; 
      std::cout << "\n" << S[i] << std::endl;
      std::cout << "\n" 
        << itpp::round_to_zero(v - S[i],thresh) << std::endl; 
      return false; 
    }
  return true; 
}

int 
main(void)
{

  // pars
  int N(6) , M(3) , B(50); 

  // j stands for junk
  itpp::Mat<std::complex<double> > jN, jM, A; 
  itpp::Mat<std::complex<double> > jUN,jVN,jUM,jVM,U,V; 
  itpp::Vec<double> jSN,jSM,S,Sinv; 

  jN = itpp::randn_c(N,N);
  itpp::svd(jN,jUN,jSN,jVN); 

  
  jM = itpp::randn_c(M,M);
  itpp::svd(jM,jUM,jSM,jVM); 

 
  // targets
  U = jUN; 
  U.del_cols(M,N-1); 
  V = jVM; 
  S = itpp::randu(M);  // (0,1)
  Sinv = S; 
  Sinv.zeros(); 
  for(int i =0; i < S.size(); ++i)
   Sinv[i] = 1./S[i];

  A = U * itpp::diag(S) * itpp::hermitian_transpose(V);  

  itpp::svd(A,jUN,jSN,jVN); 
  std::cout << "\n\n raw " << std::endl;
  std::cout << "\n" << itpp::round_to_zero(jUN,1e-6) << std::endl;
  std::cout << "\n" << itpp::round_to_zero(jSN,1e-6) << std::endl;
  std::cout << "\n" << itpp::round_to_zero(jVN,1e-6) << std::endl;



  // fake up ensembles
  std::vector<itpp::Mat<std::complex<double> > > vA(B,A), vU(B,U), vV(B,V); 
  std::vector<itpp::Vec<double> > vS(B,S); 
  


  SEMBLE::SembleMatrix<std::complex<double> > eA(vA) , eU(vU), eV(vV), cSinv; 
  SEMBLE::SembleMatrix<double> eSinv; 
  SEMBLE::SembleVector<double> eS(vS); 

  // std::string log = svdNonSquare_THREADED(eA,eU,eS,eV); 
  std::string log = svdNonSquare(eA,eU,eS,eV); 

  std::cout << log << std::endl; 

  SEMBLE::SembleMatrix<std::complex<double> > cS;
  cS = SEMBLE::recast<std::complex<double>,double>(SEMBLE::diag(eS) ); 

  // test if we can cut it up and put it back 
  test_equiv(A,eU * cS * SEMBLE::adj(eV) );

  itpp::Mat<std::complex<double> > testI; 
  SEMBLE::SembleMatrix<std::complex<double> > eTestI; 

  SEMBLE::pseudoInvertCond(SEMBLE::diag(eS),eSinv); 
  cSinv = SEMBLE::recast<std::complex<double>,double>(eSinv); 

  testI =  V * itpp::diag(Sinv) * itpp::hermitian_transpose(U) * A;
  eTestI =  eV * cSinv * SEMBLE::adj(eU) * eA ; 
    
  // test the inversion
  std::cout << "should be identity \n" 
    << itpp::round_to_zero( eTestI.mean() , 1e-6) << std::endl; 

  return 0;
}
