/*! \file stfinvfourier.cc
 * \brief a base class for all engines which operate in the Fourier domain (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * \author Thomas Forbriger
 * \date 08/05/2011
 * 
 * a base class for all engines which operate in the Fourier domain (implementation)
 * 
 * Copyright (c) 2011 by Thomas Forbriger (BFO Schiltach) 
 *
 * ----
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * ----
 *
 * 
 * REVISIONS and CHANGES 
 *  - 08/05/2011   V1.0   Thomas Forbriger
 *  - 30/09/2011   V1.1   implemented handling of additional time series pairs
 *  - 04/10/2011   V1.2   correction in debug message
 *  - 14/10/2015   V1.3   new end-user usage functions
 * 
 * ============================================================================
 */
#define STFINV_STFINVFOURIER_CC_VERSION \
  "STFINV_STFINVFOURIER_CC   V1.3"

#include <sstream>
#include <stfinv/stfinvfourier.h>
#include <stfinv/stfinvfourier_summary_usage.h>
#include <stfinv/stfinvfourier_description_usage.h>
#include <stfinv/debug.h>
#include <aff/subarray.h>
#include <aff/slice.h>
#include <aff/arrayoperators.h>

namespace stfinv {

  /*! 
   * - References to signal storage are available
   * - Parameters are parsed and available through the base class
   * What has to be done here:
   * -# The FFT engine must be initialized appropriately
   *    -# Time domain padding has to be read from the parameters
   *    -# Number of sample modulo has to be read from parameters
   *    -# number of samples must be calculated from these values
   *    -# based on these values workspace arrays must be constructed
   *    -# The FFT engine must be constructed from the workspace arrays.
   *       This requires a copy constructor for the FFT engine, which
   *       currently is not present.
   *       Alternatively, we could only use a pointer to the engine.
   */
  STFFourierDomainEngine::STFFourierDomainEngine(const stfinv::Tvectoroftriples& triples,
                                                 const stfinv::Waveform& stf,
                                                 const stfinv::Tvectorofpairs& pairs,
                                                 const std::string& parameters)
    :Tbase(triples, stf, pairs, parameters)
  { 
    this->initialize();
  } // STFFourierDomainEngine::STFFourierDomainEngine

  /*----------------------------------------------------------------------*/

  /*! 
   * - References to signal storage are available
   * - Parameters are parsed and available through the base class
   * What has to be done here:
   * -# The FFT engine must be initialized appropriately
   *    -# Time domain padding has to be read from the parameters
   *    -# Number of sample modulo has to be read from parameters
   *    -# number of samples must be calculated from these values
   *    -# based on these values workspace arrays must be constructed
   *    -# The FFT engine must be constructed from the workspace arrays.
   *       This requires a copy constructor for the FFT engine, which
   *       currently is not present.
   *       Alternatively, we could only use a pointer to the engine.
   */
  STFFourierDomainEngine::STFFourierDomainEngine(const stfinv::Tvectoroftriples& triples,
                                                 const stfinv::Waveform& stf,
                                                 const std::string& parameters)
    :Tbase(triples, stf, parameters)
  { 
    this->initialize();
  } // STFFourierDomainEngine::STFFourierDomainEngine

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::help(std::ostream& os) const
  {
    STFFourierDomainEngine::classhelp(os);
  } // void STFFourierDomainEngine::help(std::ostream& os) const

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::usage(std::ostream& os) const
  {
    STFFourierDomainEngine::classusage(os);
  } // void STFFourierDomainEngine::usage(std::ostream& os) const


  /*----------------------------------------------------------------------*/

  const char* STFFourierDomainEngine::name() const
  {
    return("STFFourierDomainEngine");
  } //  const char const* STFFourierDomainEngine::name() const

  /*----------------------------------------------------------------------*/

  /*! \brief online help text giving information on options
   *
   * This must be kept synchronous with the options used by
   * STFFourierDomainEngine::initialize()
   */
  void STFFourierDomainEngine::classhelp(std::ostream& os)
  {
    os << stfinvfourier_summary_usage;
  } // void STFFourierDomainEngine::classhelp(std::ostream& os)

  /*----------------------------------------------------------------------*/

  /*! \brief online help text giving information on options
   *
   * This must be kept synchronous with the options used by
   * STFFourierDomainEngine::initialize()
   */
  void STFFourierDomainEngine::classusage(std::ostream& os)
  {
    os << stfinvfourier_description_usage;
    os << std::endl;
    Tbase::classusage(os);
  } // void STFFourierDomainEngine::classusage(std::ostream& os)

  /*----------------------------------------------------------------------*/

  /*! \brief Create FFT engines
   *
   * \par Number of samples
   * For the FFT we usually use a number of samples larger than that of the
   * underlying time series. 
   * This way we at least partially avoid the undesired effect of the cyclic
   * discrete Fourier transform, which can lead to wrap-around.
   * Currently there are three ways to modify the number of samples actually
   * used:
   * -# Option \c fpad specifies an padding factor. 
   *    The number of samples simply is larger than the number of samples of
   *    the underlying time series by this factor.
   * -# Option \c fpow2:
   *    If set, the next power of two larger than the fpad*nsamples is used.
   * -# Option \c fdiv specifies a divisor.
   *    If set and if \c fpow2 is not set the number of samples will be the
   *    next integer multiple of the divisor larger than fpad*nsamples.
   *
   * \par Workspace
   * Two FFT engines will be created:
   * -# One engine (\c STFFourierDomainEngine::Mfftengineinput) being shared
   *  by recorded data and synthetic data, because both have to be transformed
   *  to Fourier domain at once.
   * -# One engine (\c STFFourierDomainEngine::Mfftengineoutput) being shared
   *  by the stf and the convolved synthetics, because both have to be
   *  transformed to the time domain at once.
   */
  void STFFourierDomainEngine::initialize()
  {
    // extract parameter values
    // ------------------------
    // padding factor
    double padfactor;
    {
      std::istringstream is(this->parameter("fpad","1.5"));
      is >> padfactor;
    }
    STFINV_assert(padfactor >= 1.,
                  "ERROR: parameter for option \"fpad\" not larger or equal 1");
    // flag: use power of two
    bool poweroftwo=(this->parameter("fpow2","false")=="true");
    // flag: use integer multiples of divisor
    bool divisorset=this->parameterisset("fdiv");
    // number of samples shall be integer power of divisor
    unsigned int divisor=1;
    if (divisorset)
    {
      std::istringstream is (this->parameter("fdiv","100")); 
      is >> divisor;
      STFINV_assert(divisor > 0,
                  "ERROR: parameter for option \"fdiv\" not larger than 0");
    }

    // define number of samples to be used by Fourier engine
    // -----------------------------------------------------
    // use at least the number of samples sepcified by the padfactor
    unsigned int nsamples=padfactor*this->nsamples();
    if (nsamples<this->nsamples()) { nsamples=this->nsamples(); }
    // power of two has precendence
    if (poweroftwo)
    {
      unsigned int n=2;
      while (n<nsamples) { n *= 2; }
      nsamples=n;
    }
    else if (divisorset)
    {
      unsigned int rest=nsamples % divisor;
      nsamples -= rest;
      nsamples += divisor;
    } // if (poweroftwo) or (divisorset)

    // allocate workspace by creating engines
    // --------------------------------------
    Mfftengineinput=fourier::fft::DRFFTWAFFArrayEngine(
                                     2*this->nreceivers()+this->npairs(),
                                                       nsamples);
    Mfftengineoutput=fourier::fft::DRFFTWAFFArrayEngine(
                                     1+this->nreceivers()+this->npairs(),
                                                        nsamples);

    // set time shift for STF if requested
    // -----------------------------------
    {
      std::istringstream is(this->parameter("tshift","0."));
      is >> Mtshift;
    }
    Mapplyshift=this->parameterisset("tshift");

    STFINV_debug(Mdebug&1, "STFFourierDomainEngine::initialize()",
                 "this->nsamples()=" << this->nsamples() << " " <<
                 "padfactor=" << padfactor << " " <<
                 "divisor=" << divisor << " " <<
                 "nsamples=" << nsamples << " " <<
                 "Mapplyshift=" << Mapplyshift << " " <<
                 "Mtshift=" << Mtshift << " ");
  } // void STFFourierDomainEngine::initialize() 

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::fftinput()
  {
    this->getinput();
    Mfftengineinput.r2c();
  } // void STFFourierDomainEngine::fftinput() 

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::fftoutput()
  {
    this->convolve();
    if (Mapplyshift) { this->stfshift(); }
    Mfftengineoutput.c2r();
    this->putoutput();
  } // void STFFourierDomainEngine::fftoutput() 

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum
    STFFourierDomainEngine::recordingspec() const
  {
    TAspectrum inspecarray=Mfftengineinput.spectrum();
    TAspectrum subarray=aff::subarray(inspecarray)()(0,this->nreceivers()-1);
    subarray.shape().setfirst(0,0);
    subarray.shape().setfirst(1,0);
    return(subarray);
  } // STFFourierDomainEngine::TAspectrum STFFourierDomainEngine::recordingspec() const

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum
    STFFourierDomainEngine::syntheticspec() const
  {
    TAspectrum inspecarray=Mfftengineinput.spectrum();
    TAspectrum subarray=aff::subarray(inspecarray)()(this->nreceivers(),
                                       (2*this->nreceivers())-1);
    subarray.shape().setfirst(0,0);
    subarray.shape().setfirst(1,0);
    return(subarray);
  } // STFFourierDomainEngine::TAspectrum STFFourierDomainEngine::syntheticspec() const

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum STFFourierDomainEngine::stfspec() const
  {
    return(Mfftengineoutput.spectrum(this->nreceivers()));
  } // STFFourierDomainEngine::Tspectrum STFFourierDomainEngine::stfspec() const

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum 
    STFFourierDomainEngine::recordingcoeff(const unsigned int& i) const
  {
    STFFourierDomainEngine::TAspectrum rec=this->recordingspec();
    STFFourierDomainEngine::TAspectrum slice=aff::slice(rec)(0,i);
    slice.shape().setfirst(0,0);
    slice.shape().setfirst(1,0);
    return(slice);
  } // STFFourierDomainEngine::recordingcoeff(const unsigned int& i) const

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum 
    STFFourierDomainEngine::syntheticcoeff(const unsigned int& i) const
  {
    STFINV_debug(Mdebug&8, "STFFourierDomainEngine::syntheticcoeff",
                 "i=" << i);
    STFFourierDomainEngine::TAspectrum syn=this->syntheticspec();
    STFFourierDomainEngine::TAspectrum slice=aff::slice(syn)(0,i);
    slice.shape().setfirst(0,0);
    slice.shape().setfirst(1,0);
    STFINV_debug(Mdebug&8, "STFFourierDomainEngine::syntheticcoeff",
                 "slice.first(0)=" << slice.first(0) << " " <<
                 "slice.last(0)=" << slice.last(0) << " "
                 "slice.first(1)=" << slice.first(1) << " " <<
                 "slice.last(1)=" << slice.last(1) << " ");
    return(slice);
  } // STFFourierDomainEngine::syntheticcoeff(const unsigned int& i) const

  /*----------------------------------------------------------------------*/

  STFFourierDomainEngine::TAspectrum::Tvalue& 
    STFFourierDomainEngine::stfcoeff(const unsigned int& i) const
  {
    STFFourierDomainEngine::TAspectrum stffft=this->stfspec();
    return(stffft(i));
  } // STFFourierDomainEngine::stfcoeff(const unsigned int& i) const

  /*----------------------------------------------------------------------*/

  double STFFourierDomainEngine::frequency(const unsigned int& i) const
  {
    return(static_cast<double>(i)/(this->dt()*Mfftengineinput.nsamples()));
  } // double STFFourierDomainEngine::frequency(const unsigned int& i) const

  /*----------------------------------------------------------------------*/

  unsigned int STFFourierDomainEngine::nfreq() const
  {
    return(Mfftengineinput.nfrequencies());
  } // unsigned int STFFourierDomainEngine::nfreq() const

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::getinput()
  {
    // clear workspace through reference
    TAseries sarray=Mfftengineinput.series();
    sarray=0.;

    // cycle through receivers
    for (unsigned int i=0; i<this->nreceivers(); ++i)
    {
      // get references to time series in workspace
      TAseries recordingref=Mfftengineinput.series(i);
      TAseries syntheticref=Mfftengineinput.series(i+this->nreceivers());
      // copyin function copies as many elements as possible
      recordingref.copyin(this->recording(i));
      syntheticref.copyin(this->synthetic(i));
    } // for (unsigned int i=0; i<this->nreceivers(); ++i)

    // cycle through additional time series pairs
    for (unsigned int i=0; i<this->npairs(); ++i)
    {
      // get references to time series in workspace
      TAseries syntheticref=Mfftengineinput.series(i+2*this->nreceivers());
      // copyin function copies as many elements as possible
      syntheticref.copyin(this->series(i));
    } // for (unsigned int i=0; i<this->npairs(); ++i)
  } // void STFFourierDomainEngine::getinput() 

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::putoutput()
  {
    // scaling factor for source correction filter
    double scalingstf=Mfftengineoutput.scale_series(this->dt());
    // scaling factor for convolved synthetics
    double scalingcs=Mfftengineoutput.scale_spectrum(this->dt())
      *scalingstf;

    // cycle through receivers
    for (unsigned int i=0; i<this->nreceivers(); ++i)
    {
      // get references to time series 
      stfinv::Tseries convolvedsynthetics=this->convolvedsynthetic(i);
      // copyin function copies as many elements as possible
      convolvedsynthetics.copyin(Mfftengineoutput.series(i)*scalingcs);
    } // for (unsigned int i=0; i<this->nreceivers(); ++i)
      
    // copy stf too; get reference to series and use copyin
    stfinv::Tseries stf=this->stf();
    stf.copyin(Mfftengineoutput.series(this->nreceivers())
               *scalingstf);

    // cycle through additional time series pairs
    for (unsigned int i=0; i<this->npairs(); ++i)
    {
      // get references to time series 
      stfinv::Tseries convolvedseries=this->convolvedseries(i);
      // copyin function copies as many elements as possible
      convolvedseries.copyin(Mfftengineoutput.series(1+i+this->nreceivers())
                             *scalingcs);
    } // for (unsigned int i=0; i<this->npairs(); ++i)

  } // void STFFourierDomainEngine::putoutput() 

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::convolve()
  {
    TAspectrum synthetic=this->syntheticspec();
    TAspectrum convolvedsynthetic
      =aff::subarray(Mfftengineoutput.spectrum())()(0,this->nreceivers()-1);
    convolvedsynthetic.shape().setfirst(0,0);
    convolvedsynthetic.shape().setfirst(1,0);
    TAspectrum stfcoeff=this->stfspec();

    // cycle through all receivers
    for (unsigned int i=0; i<this->nreceivers(); ++i)
    {
      for (unsigned int j=0; j<this->nfreq(); ++j)
      {
        convolvedsynthetic(j,i)=synthetic(j,i)*stfcoeff(j);
      } // for (unsigned int j=0; j<this->nfreq(); ++j)
    } // for (unsigned int i=0; i<this->nreceivers(); ++i)

    // copy signals for additional pairs
    if (this->npairs()>0)
    {
      // reference to input series Fourier coefficients of additional pairs
      TAspectrum inspecarray=Mfftengineinput.spectrum();
      TAspectrum seriesspecarray
        =aff::subarray(inspecarray)()(2*this->nreceivers(),
                                         (2*this->nreceivers()
                                          +this->npairs()-1));
      seriesspecarray.shape().setfirst(0,0);
      seriesspecarray.shape().setfirst(1,0);

      // reference to convolved series Fourier coefficients of additional
      // pairs
      TAspectrum convolvedseries
        =aff::subarray(Mfftengineoutput.spectrum())()(this->nreceivers()+1,
                                                      this->nreceivers()
                                                      +this->npairs());
      convolvedseries.shape().setfirst(0,0);
      convolvedseries.shape().setfirst(1,0);
      // cycle through all pairs
      for (unsigned int i=0; i<this->npairs(); ++i)
      {
        for (unsigned int j=0; j<this->nfreq(); ++j)
        {
          convolvedseries(j,i)=seriesspecarray(j,i)*stfcoeff(j);
        } // for (unsigned int j=0; j<this->nfreq(); ++j)
      } // for (unsigned int i=0; i<this->npairs(); ++i)
    } // if (this->npairs()>0)

  } // void STFFourierDomainEngine::convolve() 

  /*----------------------------------------------------------------------*/

  void STFFourierDomainEngine::stfshift()
  {
    STFINV_debug(Mdebug&2, "STFFourierDomainEngine::stfshift()",
                 "apply time shift of " << Mtshift << " seconds");
    const STFFourierDomainEngine::TAspectrum::Tvalue sfac
      =-STFFourierDomainEngine::TAspectrum::Tvalue(0.,1.)
      *3.141592653589793115998*2.*Mtshift;
    // STFINV_DUMP( sfac );
    STFFourierDomainEngine::TAspectrum stffft=this->stfspec();
    for (unsigned int j=0; j<this->nfreq(); ++j)
    {
      stffft(j) *= std::exp(sfac*this->frequency(j));
      // STFINV_DUMP( std::exp(sfac*this->frequency(j)) );      
    } // for (unsigned int j=0; j<this->nfreq(); ++j)
  } // void STFFourierDomainEngine::stfshift() 

} // namespace stfinv

/* ----- END OF stfinvfourier.cc ----- */
