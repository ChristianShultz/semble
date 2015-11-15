#ifndef SEMBLE_ALGEBRA_H_H_GUARD
#define SEMBLE_ALGEBRA_H_H_GUARD

#include "semble_matrix.h"
#include "semble_vector.h"
#include "semble_meta.h"


namespace SEMBLE
{
  //Semble Vector Algebra
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //rescale
  template<class T>
    SembleVector<T> rescaleEnsemDown(const SembleVector<T> &V)
    {
      SembleVector<T> dum(V);
      dum.rescaleEnsemDown();
      return dum;
    }

  //rescale
  template<class T>
    SembleVector<T> rescaleEnsemUp(const SembleVector<T> &V)
    {
      SembleVector<T> dum(V);
      dum.rescaleEnsemUp();
      return dum;
    }

  //equivilance
  template<class T>
    bool operator==(const SembleVector<T> &V1, const SembleVector<T> &V2)
    {
      if(V1.getB() != V2.getB())
      {
        return false;
      }
      else
      {
        int dim = V1.getB();

        for(int bin = 0; bin < dim; ++bin)
        {
          if(V1[bin] != V2[bin])
            return false;
        }
      }

      return true;
    }

  //negation equivilance
  template<class T>
    bool operator!=(const SembleVector<T> &V1, const SembleVector<T> &V2)
    {
      return !(V1 == V2);
    }

  //addition
  template<class T>
    SembleVector<T> operator+(const SembleVector<T> &V1, const SembleVector<T> &V2)
    {
      SembleVector<T> Vr(V1);
      Vr += V2;
      return Vr;
    }

  template<class T>
    SembleVector<T> operator+(const SembleVector<T> &V1, const itpp::Vec<T> &V2)
    {
      SembleVector<T> dum(V1);
      dum += V2;
      return dum;
    }

  template<class T>
    SembleVector<T> operator+(const itpp::Vec<T> &V1, const SembleVector<T> &V2)
    {
      return V2 + V1;
    }

  //subtraction
  template<class T>
    SembleVector<T> operator-(const SembleVector<T> &V1, const SembleVector<T> &V2)
    {
      SembleVector<T> dum(V1);
      dum -= V2;
      return dum;
    }

  template<class T>
    SembleVector<T> operator-(const SembleVector<T> &V1, const itpp::Vec<T> &V2)
    {
      return V1 + (-1.0 * V2);
    }

  template<class T>
    SembleVector<T> operator-(const itpp::Vec<T> &V1, const SembleVector<T> &V2)
    {
      return (-V2) + V1;
    }

  //multiplication
  template<class T>
    SembleVector<T> operator*(const T &lhs, const SembleVector<T> &V)
    {
      SembleVector<T> dum(V);
      dum *= lhs;
      return dum;
    }

  template<class T>
    SembleVector<T> operator*(const SembleVector<T> &V, const T &rhs)
    {
      return rhs * V;
    }

  template<class T>
    SembleVector<T> operator*(const SembleVector<T> &V, const typename PromoteScalar<T>::Type &rhs)
    {
      SembleVector<T> dum(V);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleVector<T> operator*(const typename PromoteScalar<T>::Type &lhs, const SembleVector<T> &V)
    {
      return V * lhs;
    }

  template<class T>
    SembleVector<T> operator*(const SembleVector<T> &V, const typename PromoteEnsem<T>::Type &rhs)
    {
      SembleVector<T> dum(V);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleVector<T> operator*(const typename PromoteEnsem<T>::Type &lhs, const SembleVector<T> &V)
    {
      return V * lhs;
    }

  //dot product -- NB this does the CC
  template<class T>  
    typename PromoteEnsem<T>::Type operator*(const SembleVector<T> &lhs, const SembleVector<T> &rhs)
    {
      SembleVector<T> dum(rhs);
      dum.conj();
      typename PromoteEnsem<T>::Type ddum = dum.dot(lhs);
      return ddum;
    }

  template<class T>
    SembleVector<T> operator*(const itpp::Mat<T> &M, const SembleVector<T> &rhs)
    {
      SembleVector<T> dum(rhs);
      int nb = dum.getB();
      for(int v = 0; v < nb; ++v)
        dum[v] = M*rhs[v];
      return dum;
    }

  template<class T>
    SembleVector<T> operator*(const SembleVector<T> &lhs, const itpp::Mat<T> &M)
    {
      SembleVector<T> dum(lhs);
      dum.conj();
      itpp::Mat<T> Md(M.H());

      return (Md*dum).conj();
    }

  //division
  template<class T>
    SembleVector<T> operator/(const SembleVector<T> &lhs, const T &rhs)
    {
      SembleVector<T> dum(lhs);
      return dum /= rhs;
    }

  template<class T>
    SembleVector<T> operator/(const SembleVector<T> &lhs, const typename PromoteScalar<T>::Type &rhs)
    {
      SembleVector<T> dum(lhs);
      return dum /= rhs;
    }

  template<class T>
    SembleVector<T> operator/(const SembleVector<T> &lhs, const typename PromoteEnsem<T>::Type &rhs)
    {
      SembleVector<T> dum(lhs);
      return dum /= rhs;
    }

  template<class T>
    itpp::Vec<T> mean(const SembleVector<T> &in)
    {
      return in.mean();
    }

  template<class T>
    itpp::Vec<T> variance(const SembleVector<T> &in)
    {
      return in.variance();
    }

  template<class T>
    SembleVector<T> complexConjugate(const SembleVector<T> &in)
    {
      SembleVector<T> dum(in.getB(), in.getN());
      int dim = dum.getB();

      for(int bin = 0; bin < dim; ++bin)
        dum[bin] = itpp::conj(in[bin]);

      return dum;
    }

  //Semble Matrix Algebra
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //rescale
  template<class T>
    SembleMatrix<T> rescaleEnsemDown(const SembleMatrix<T> &M)
    {
      SembleMatrix<T> dum(M);
      dum.rescaleSembleDown();
      return dum;
    }

  template<class T>
    SembleMatrix<T> rescaleEnsemUp(const SembleMatrix<T> &M)
    {
      SembleMatrix<T> dum(M);
      dum.rescaleSembleUp();
      return dum;
    }

  //equivilance
  template<class T>
    bool operator==(const SembleMatrix<T> &M1, const SembleMatrix<T> &M2)
    {
      if(M1.getB() != M2.getB())
      {
        return false;
      }
      else
      {
        int dim = M1.getB();

        for(int bin = 0; bin < dim; ++bin)
          if(M1[bin] != M2[bin])
            return false;
      }

      return true;
    }

  //negation of equivilance
  template<class T>
    bool operator!=(const SembleMatrix<T> &lhs, const SembleMatrix<T> &rhs)
    {
      return !(lhs == rhs);
    }

  //addition
  template<class T>
    SembleMatrix<T> operator+(const SembleMatrix<T> &lhs, const SembleMatrix<T> &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum += rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator+(const SembleMatrix<T> &lhs, const itpp::Mat<T> &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum += rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator+(const itpp::Mat<T> &lhs, const SembleMatrix<T> &rhs)
    {
      return (rhs + lhs);
    }

  //subtraction
  template<class T>
    SembleMatrix<T> operator-(const SembleMatrix<T> &lhs, const SembleMatrix<T> &rhs)
    {
      return (lhs + (-rhs));
    }

  template<class T>
    SembleMatrix<T> operator-(const SembleMatrix<T> &lhs, const itpp::Mat<T> &rhs)
    {
      return (lhs + (-rhs));
    }

  template<class T>
    SembleMatrix<T> operator-(const itpp::Mat<T> &lhs, const SembleMatrix<T> &rhs)
    {
      return (-(rhs + (-lhs)));
    }

  //multiplication
  template<class T>
    SembleMatrix<T> operator*(const SembleMatrix<T> &lhs, const T &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator*(const T &lhs, const SembleMatrix<T> &rhs)
    {
      return rhs * lhs;
    }

  template<class T>
    SembleMatrix<T> operator*(const SembleMatrix<T> &lhs, const typename PromoteScalar<T>::Type &rhs)
    {
      return lhs * toScalar(rhs);
    }

  template<class T>
    SembleMatrix<T> operator*(const typename PromoteScalar<T>::Type &lhs, const SembleMatrix<T> &rhs)
    {
      return rhs * toScalar(lhs);
    }

  template<class T>
    SembleMatrix<T> operator*(const SembleMatrix<T> &lhs, const typename PromoteEnsem<T>::Type &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator*(const typename PromoteEnsem<T>::Type &lhs, const SembleMatrix<T> &rhs)
    {
      return rhs * lhs;
    }

  template<class T>
    SembleMatrix<T> operator*(const SembleMatrix<T> &lhs, const SembleMatrix<T> &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator*(const SembleMatrix<T> &lhs, const itpp::Mat<T> &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum *= rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator*(const itpp::Mat<T> &lhs, const SembleMatrix<T> &rhs)
    {
      SembleMatrix<T> dum(adj(rhs));
      dum *= itpp::hermitian_transpose(lhs);
      return adj(dum);
    }

  //division
  template<class T>
    SembleMatrix<T> operator/(const SembleMatrix<T> &lhs, const T &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum /= rhs;
      return dum;
    }

  template<class T>
    SembleMatrix<T> operator/(const SembleMatrix<T> &lhs, const typename PromoteScalar<T>::Type &rhs)
    {
      return (lhs / toScalar(rhs));
    }

  template<class T>
    SembleMatrix<T> operator/(const SembleMatrix<T> &lhs, const typename PromoteEnsem<T>::Type &rhs)
    {
      SembleMatrix<T> dum(lhs);
      dum /= rhs;
      return dum;
    }

  template<class T>
    itpp::Mat<T> mean(const SembleMatrix<T> &in)
    {
      return in.mean();
    }

  template<class T>
    itpp::Mat<T> variance(const SembleMatrix<T> &in)
    {
      return in.variance();
    }

  template<class T>
    SembleMatrix<T> complexConjugate(const SembleMatrix<T> &in)
    {
      SembleMatrix<T> dum(in.getB(), in.getN(), in.getM());
      int dim = dum.getB();

      for(int bin = 0; bin < dim; ++bin)
        dum[bin] = itpp::conj(in[bin]);

      return dum;
    }


  // an asymmetric cast if underlying type is possible
  // this will only work if the toScalar func is supported for the underlying in 
  // type and there is an implicit cast available from in_type to out_type
  template<typename outT, typename inT>
    SEMBLE::SembleMatrix<outT> recast(const SEMBLE::SembleMatrix<inT> &i)
    {
      const int B = i.getB();
      const int N = i.getN();
      const int M = i.getM();

      SEMBLE::SembleMatrix<outT> o(B,N,M);
      for(int b = 0; b < B; ++b)
        for(int n = 0; n < N; ++n)
          for(int m = 0; m < M; ++m)
            o.setElement(b,n,m,outT(SEMBLE::toScalar(i(b,n,m))));
      return o;
    }


  SEMBLE::SembleVector<std::complex<double> > round_to_zero(const SEMBLE::SembleVector<std::complex<double> > &in, const double thresh=1e-14); 
  SEMBLE::SembleVector<double> round_to_zero(const SEMBLE::SembleVector<double> &in, const double thresh=1e-14);

}
#endif
