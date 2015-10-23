/*! \file stfinvfourier.h
 * \brief a base class for all engines which operate in the Fourier domain (prototypes)
 * 
 * ----------------------------------------------------------------------------
 * 
 * \author Thomas Forbriger
 * \date 08/05/2011
 * 
 * a base class for all engines which operate in the Fourier domain (prototypes)
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
 * REVISIONS and CHANGES 
 *  - 08/05/2011   V1.0   Thomas Forbriger
 *  - 30/09/2011   V1.1   implemented handling of additional time series pairs
 *  - 14/10/2015   V1.2   new end-user usage functions
 * 
 * ============================================================================
 */

// include guard
#ifndef STFINV_STFINVFOURIER_H_VERSION

#define STFINV_STFINVFOURIER_H_VERSION \
  "STFINV_STFINVFOURIER_H   V1.2"

#include <stfinv/stfinvbase.h>
#include <aff/array.h>
#include <fourier/fftwaffar.h>

namespace stfinv {

  /*! \brief Base class for all engines which operate in the Fourier domain
   *
   * This is just a base class.
   * The constructor is protected and should only be called from a derived
   * class.
   *
   * This class maintains a workspace for Fourier transforms.
   * It provides the FFT from input signals to the workspace through a member
   * functions as well as the convolution of the synthetic data with a given
   * source wavelet Fourier transform and a subsequent FFT to time domain for
   * the convolved synthetics as well as the source correction filter separately.
   *
   * \par What STFFourierDomainEngine does for you
   * All derived classes call STFFourierDomainEngine::fftinput prior to
   * processing.
   * This function copies (by calling STFFourierDomainEngine::getinput) the
   * users workspace to the Fourier transform engine
   * STFFourierDomainEngine::Mfftengineinput and transforms the time series to
   * the Fourier domain.
   * 
   * \par
   * When processing has finished, the derived classes should call
   * STFFourierDomainEngine::fftoutput.
   * This function convolves the synthetic data with the source correction
   * filter (STFFourierDomainEngine::convolve).
   * Then it applies a time shift to the source correction filter if requested
   * (STFFourierDomainEngine::stfshift).
   * The convolved synthetics as well as the source correction filter then are
   * transformed to time domain and written to the users workspace
   * ((STFFourierDomainEngine::putoutput).
   *
   * \par Layout of Fourier transform arrays
   * The workspace for the Fourier transform engine is initialized by
   * fourier::fft::DRFFTWAFFArrayEngine.
   */
  class STFFourierDomainEngine: public stfinv::STFBaseEngine {
    public:
      //! \brief typedef to refer to base class
      typedef stfinv::STFBaseEngine Tbase;
      //! \brief type of underlying Fourier engine
      typedef fourier::fft::DRFFTWAFFArrayEngine Tfftengine;
      //! \brief type of array for time series values
      typedef Tfftengine::TAseries TAseries;
      //! \brief type of array for Fourier transforms
      typedef Tfftengine::TAspectrum TAspectrum;
    protected:
      /*! \brief Constructor.
       */
      STFFourierDomainEngine(const stfinv::Tvectoroftriples& triples,
                             const stfinv::Waveform& stf,
                             const std::string& parameters);
      /*! \brief Constructor.
       */
      STFFourierDomainEngine(const stfinv::Tvectoroftriples& triples,
                             const stfinv::Waveform& stf,
                             const stfinv::Tvectorofpairs& pairs,
                             const std::string& parameters);
    public:
      //! \brief abstract base requires virtual destructor
      virtual ~STFFourierDomainEngine() { }
      //! \brief Start engine 
      virtual void exec() { STFINV_baseillegal; }
      //! \brief print online help
      virtual void help(std::ostream& os=std::cout) const;
      //! \brief print online help
      static void classhelp(std::ostream& os=std::cout);
      //! \brief print detailed description
      virtual void usage(std::ostream& os=std::cout) const;
      //! \brief print detailed description
      static void classusage(std::ostream& os=std::cout);
      //! \brief return name of engine
      virtual const char* name() const;
    protected:
      /*! \brief copy input signals to workspace and
       * transform input workspace to Fourier domain
       */
      void fftinput();
      /*! \brief convolve synthetics with Fourier transform of stf and
       * transform convolved synthetics and stf to time domain and pass
       * signals to user memory space
       */
      void fftoutput();
      //! \brief return reference to Fourier transform of recorded data
      TAspectrum recordingspec() const;
      //! \brief return reference to Fourier transform of synthetics
      TAspectrum syntheticspec() const;
      //! \brief return reference to Fourier transform of stf
      TAspectrum stfspec() const;
      /*! \brief return reference to Fourier coeffients of recorded data for
       * frequency i
       */
      TAspectrum recordingcoeff(const unsigned int& i) const;
      /*! \brief return reference to Fourier coefficients of synthetics for
       * frequency i
       */
      TAspectrum syntheticcoeff(const unsigned int& i) const;
      /*! \brief return reference to Fourier coefficients of stf for frequency
       * i
       */
      TAspectrum::Tvalue& stfcoeff(const unsigned int& i) const;
      //! \brief return value of frequency i in Hz
      double frequency(const unsigned int& i) const;
      //! \brief return number of frequencies in use
      unsigned int nfreq() const;
    private:
      //! \brief initialize work space
      void initialize();
      /*! \brief copy input time series for recorded data and synthetics
       * to workspace
       */
      void getinput();
      /*! \brief copy workspace time series for convolved synthetics and stf
       * to user memory
       */
      void putoutput();
      /*! \brief convolve synthetics with stf
       */
      void convolve();
      /*! \brief apply time shift to stf prior to FFT to time domain
       */
      void stfshift();
      
    // member data
    // -----------
    private:
      /*! \brief combined FFT engine for recorded data and synthetics
       * and additional time series 
       *
       * Layout of workspace:
       * <TABLE>
       * <TR>
       *   <TH> series index </TH>
       *   <TH> type of series </TH>
       * </TR>
       * <TR>
       *   <TD> 0 ... N-1 </TD>
       *   <TD> recorded data for N receivers </TD>
       * </TR>
       * <TR>
       *   <TD> N ... 2*N-1 </TD>
       *   <TD> synthetic data for N receivers </TD>
       * </TR>
       * <TR>
       *   <TD> 2*N ... 2*N+M-1 </TD>
       *   <TD> input series for M pairs </TD>
       * </TR>
       * </TABLE>
       * N is returned by function nreceivers().
       * M is returned by function npairs().
       */
      Tfftengine Mfftengineinput;
      /*! \brief combined FFT engine for stf and convolved synthetics
       * and additional convolved time series 
       *
       * Layout of workspace:
       * <TABLE>
       * <TR>
       *   <TH> series index </TH>
       *   <TH> type of series </TH>
       * </TR>
       * <TR>
       *   <TD> 0 ... N-1 </TD>
       *   <TD> convolved synthetics for N receivers </TD>
       * </TR>
       * <TR>
       *   <TD> N </TD>
       *   <TD> source correction filter </TD>
       * </TR>
       * <TR>
       *   <TD> N+1 ... N+M </TD>
       *   <TD> convolved series for M pairs </TD>
       * </TR>
       * </TABLE>
       * N is returned by function nreceivers().
       * M is returned by function npairs().
       */
      Tfftengine Mfftengineoutput;
      /*! \brief time shift to be applied to STF in order to expose
       * acausal parts
       */
      double Mtshift;
      //! \brief true if shift must be applied
      bool Mapplyshift;
  }; // class STFFourierDomainEngine

}

#endif // STFINV_STFINVFOURIER_H_VERSION (includeguard)

/* ----- END OF stfinvfourier.h ----- */
