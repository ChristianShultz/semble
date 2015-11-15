/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : semble_algebra.cc

* Purpose :

* Creation Date : 07-03-2013

* Last Modified : Thu Mar  7 16:39:19 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "semble_algebra.h"
#include "itpp/itbase.h"
#include <complex>


namespace SEMBLE
{

  SEMBLE::SembleVector<std::complex<double> > round_to_zero(const SEMBLE::SembleVector<std::complex<double> > &in, const double thresh)
  {
    SEMBLE::SembleVector<std::complex<double> > out(in); 
    itpp::Vec<std::complex<double> > mean = itpp::round_to_zero(in.mean(),thresh);
    for(int elem = 0; elem < in.getN(); ++elem)
      if(mean(elem) == std::complex<double>(0.,0.))
        out.loadEnsemElement(elem,SEMBLE::toScalar(std::complex<double>(0.,0.))*in.getEnsemElement(elem));

    return out;
  }

  SEMBLE::SembleVector<double> round_to_zero(const SEMBLE::SembleVector<double> &in, const double thresh)
  {
    SEMBLE::SembleVector<double> out(in); 
    itpp::Vec<double> mean = itpp::round_to_zero(in.mean(),thresh);
    for(int elem = 0; elem < in.getN(); ++elem)
      if(mean(elem) == 0.)
        out.loadEnsemElement(elem,SEMBLE::toScalar(0.)*in.getEnsemElement(elem));

    return out;

  }


}
