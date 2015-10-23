/*! \file stfinvany.cc
 * \brief a wrapper to any STF engine in the library (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: stfinvany.cc 4968 2013-02-01 13:58:05Z lrehor $
 * \author Thomas Forbriger
 * \date 06/05/2011
 * 
 * a wrapper to any STF engine in the library (implementation)
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
 *  - 06/05/2011   V1.0   Thomas Forbriger
 *  - 30/09/2011   V1.1   implemented handling of additional time series pairs
 *  - 04/10/2011   V1.2   renamed Fourier domain least squares engine
 * 
 * ============================================================================
 */
#define STFINV_STFINVANY_CC_VERSION \
  "STFINV_STFINVANY_CC   V1.1"
#define STFINV_STFINVANY_CC_CVSID \
  "$Id: stfinvany.cc 4968 2013-02-01 13:58:05Z lrehor $"

#include <stfinv/stfinvany.h>
#include <stfinv/stfinvfdleastsquares.h>
#include <stfinv/stfinvfixedstf.h>
#include <stfinv/stfinvidentity.h>
#include <stfinv/parameterhandler.h>
#include <stfinv/error.h>

namespace stfinv {
  
  //! \brief Constructor.
  STFEngine::STFEngine(const stfinv::Tvectoroftriples& triples,
                       const stfinv::Waveform& stf,
                       const std::string& parameters)
  {
    stfinv::Tvectorofpairs pairs;
    pairs.clear();
    this->initialize(triples, stf, pairs, parameters);
  }

  /*----------------------------------------------------------------------*/
  //! \brief Constructor.
  STFEngine::STFEngine(const stfinv::Tvectoroftriples& triples,
                       const stfinv::Waveform& stf,
                       const stfinv::Tvectorofpairs& pairs,
                       const std::string& parameters)
  {
    this->initialize(triples, stf, pairs, parameters);
  }

  /*----------------------------------------------------------------------*/

  void STFEngine::initialize(const stfinv::Tvectoroftriples& triples,
                             const stfinv::Waveform& stf,
                             const stfinv::Tvectorofpairs& pairs,
                             const std::string& parameters)
  {
    std::string para=parameters;
    std::string id=stfinv::tools::clipstring(para, ":");
    if (id == std::string(stfinv::STFEngineIdentity::ID))
    {
      STFINV_assert(pairs.size()==0,
              "ERROR: engine does not support additional time series pairs");
      Mengine=new stfinv::STFEngineIdentity(triples, stf, para);
    }
    else if (id == std::string(stfinv::STFEngineFixedWavelet::ID))
    {
      STFINV_assert(pairs.size()==0,
              "ERROR: engine does not support additional time series pairs");
      Mengine=new stfinv::STFEngineFixedWavelet(triples, stf, para);
    }
    else if ((id == std::string(stfinv::STFEngineFDLeastSquares::ID))
             || (id == std::string("fbd")))
    {
      STFINV_report_assert((id == std::string(stfinv::STFEngineFDLeastSquares::ID)),
                           "The ID \"fbd\" for this engine is deprecated",
                           "The correct ID of the Fourier domain least "
                           "squares engine is \"" <<
                           stfinv::STFEngineFDLeastSquares::ID <<
                           "\".\n"
                           "The former ID \"fbd\" may vanish in the future "
                           "and should no longer be used.");
      if (pairs.size()>0)
      {
        Mengine=new stfinv::STFEngineFDLeastSquares(triples, stf, 
                                                        pairs, para);
      }
      else
      {
        Mengine=new stfinv::STFEngineFDLeastSquares(triples, stf, para);
      }
    }
    else
    {
      std::cerr << "ERROR: engine ID " << id << " is unkown!" << std::endl;
      STFINV_abort("aborting since engine ID is not recognized");
    }
    STFINV_assert(Mengine!=0, "engine was not created correctly");
  }

  /*----------------------------------------------------------------------*/

  STFEngine::~STFEngine() 
  {
    delete Mengine;
  } // STFEngine::~STFEngine()

  /*----------------------------------------------------------------------*/

  void STFEngine::help(std::ostream& os)
  {
    os << "Currently the following engines are available:" << std::endl;
    os << "---------------------------------------------" << std::endl;
    os << std::endl;
    const int width1=10;
    /*----------------------------------------------------------------------*/
    os.width(width1); 
    os << STFEngineIdentity::ID;
    os.width(0);
    os << ": " << STFEngineIdentity::description << std::endl;
    /*----------------------------------------------------------------------*/
    os.width(width1); 
    os << STFEngineFixedWavelet::ID;
    os.width(0);
    os << ": " << STFEngineFixedWavelet::description << std::endl;
    /*----------------------------------------------------------------------*/
    os.width(width1); 
    os << STFEngineFDLeastSquares::ID;
    os.width(0);
    os << ": " << STFEngineFDLeastSquares::description << std::endl;
  } // void STFEngine::help(std::ostream& os=std::cout)

  /*======================================================================*/

  void help(std::ostream& os)
  {
    os << "Details descriptions of available engines:" << std::endl;
    os << "------------------------------------------" << std::endl;
    os << std::endl;
    STFEngineIdentity::classhelp(os);
    os << std::endl;
    STFEngineFixedWavelet::classhelp(os);
    os << std::endl;
    STFEngineFDLeastSquares::classhelp(os);

    // parameter strings
    os << "\n";
    os << "How to construct parameter strings" << std::endl;
    os << "----------------------------------" << std::endl;
    os << "A specific engine is selected by passing a parameter string.\n"
       << "This parameter string may further contain parameters to control\n"
       << "the execution mode of the engine.\n";

    os << "The parameter string starts with an ID-sequence identifying\n"
       << "the desired engine. See the list below for available engines.\n"
       << "In the parameter string the ID-sequence is terminated by a\n"
       << "colon (:).\n";

    os << "After selecting the desired engine, the interface function\n"
       << "strips of the ID-sequence as well as the colon from the\n"
       << "parameter string and initializes the engine and passes the\n"
       << "remainder of the parameter string to the engine. This\n"
       << "remainder may consist of several control parameters being\n"
       << "separated by colons (:). Each control parameter may just be\n"
       << "a flag (switch to turn an option on) or may come along with\n"
       << "a parameter value. The value of the parameter is separated\n"
       << "by an equal sign (=).\n";

    os << "\n";
    os << "Examples:\n";
    os << "- To select Fourier domain least squares and shift\n"
       << "  the returned source correction filter wavelet by 0.4s and\n"
       << "  switch on verbose mode, pass the following parameter string:\n";
    os << "    fdlsq:tshift=0.4:verbose\n";
    os << "- To select the identity engine and to switch on debug level 4:\n";
    os << "    ident:DEBUG=4\n";
    os << "- To select Fourier domain least squares, apply offset\n"
       << "  dependent weights and use a power of two to speed up the FFT:\n";
    os << "    fdlsq:pow2:exp=1.4\n";

    os << std::endl;

    STFEngine::help(os);
  } // void help(std::ostream& os=std::cout)

} // namespace stfinv

/* ----- END OF stfinvany.cc ----- */
