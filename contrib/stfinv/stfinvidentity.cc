/*! \file stfinvidentity.cc
 * \brief just find a scaling factor (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: stfinvidentity.cc 3986 2011-05-29 16:01:09Z tforb $
 * \author Thomas Forbriger
 * \date 07/05/2011
 * 
 * just find a scaling factor (implementation)
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
 *  - 07/05/2011   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define STFINV_STFINVIDENTITY_CC_VERSION \
  "STFINV_STFINVIDENTITY_CC   V1.0   "
#define STFINV_STFINVIDENTITY_CC_CVSID \
  "$Id: stfinvidentity.cc 3986 2011-05-29 16:01:09Z tforb $"

#include <stfinv/stfinvidentity.h>

namespace stfinv {

  const char* const STFEngineIdentity::ID="ident";

  const char* const STFEngineIdentity::description
    ="identity: just apply a scalar factor";

  /*----------------------------------------------------------------------*/

  void STFEngineIdentity::help(std::ostream& os) const
  {
    STFEngineIdentity::classhelp(os);
  } // void STFEngineIdentity::help(std::ostream& os) const

  /*----------------------------------------------------------------------*/

  const char* STFEngineIdentity::name() const
  {
    return("STFEngineIdentity");
  } //  const char const* STFEngineIdentity::name() const

  /*----------------------------------------------------------------------*/

  void STFEngineIdentity::initialize() 
  {
    // scale energy
    Mscaleenergy=(this->parameter("scaleenergy","false")=="true");
  } // void STFEngineIdentity::initialize()

  /*----------------------------------------------------------------------*/

  void STFEngineIdentity::exec() 
  {
    STFINV_assert(!Mscaleenergy,
                  "ERROR: energy scaling not yet implemented");
    Tseries stf=this->stf();
    stf=0.;
    stf(0)=1.;
    for (unsigned int i=0; i<this->nreceivers(); ++i)
    {
      Tseries::Tcoc synthetic=this->synthetic(i);
      Tseries convolvedsynthetic=this->convolvedsynthetic(i);
      convolvedsynthetic.copyin(synthetic);
    }
  } // void STFEngineIdentity::exec()

  /*----------------------------------------------------------------------*/

  void STFEngineIdentity::classhelp(std::ostream& os)
  {
    os << "class STFEngineIdentity (" 
      << STFEngineIdentity::ID << ")\n";
    os << STFEngineIdentity::description << "\n" << std::endl;
    os << "This engine convolves the synthetic data with a discrete delta\n"
      << "pulse so to speak. Optionally the delta-peak ist scale such that\n"
      << "the convolved synthetics will be of equal scaled energy as the\n"
      << "recordings.\n";
    os << "Options and parameters in common for Fourier engines:\n"
      << "scaleenergy   if flag is set: scale energy"
      << std::endl;
    Tbase::classhelp(os);
  } // void STFEngineIdentity::classhelp(std::ostream& os)

} // namespace stfinv

/* ----- END OF stfinvidentity.cc ----- */
