/*! \file stfinvnormalize.cc
 * \brief a Fourier domain engine which finds a normalizing source wavelet (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: stfinvnormalize.cc 4164 2011-10-04 09:00:57Z tforb $
 * \author Thomas Forbriger
 * \date 08/05/2011
 * 
 * a Fourier domain engine which finds a normalizing source wavelet (implementation)
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
 * 
 * ============================================================================
 */
#define STFINV_STFINVNORMALIZE_CC_VERSION \
  "STFINV_STFINVNORMALIZE_CC   V1.0"
#define STFINV_STFINVNORMALIZE_CC_CVSID \
  "$Id: stfinvnormalize.cc 4164 2011-10-04 09:00:57Z tforb $"

#include <stfinv/stfinvnormalize.h>

namespace stfinv {

  const char* const STFEngineNormalize::ID="fdnorm";

  const char* const STFEngineNormalize::description
    ="Fourier domain normalization";

  /*----------------------------------------------------------------------*/

  void STFEngineNormalize::help(std::ostream& os) const
  {
    STFEngineNormalize::classhelp(os);
  } // void STFEngineNormalize::help(std::ostream& os) const

  /*----------------------------------------------------------------------*/

  const char* STFEngineNormalize::name() const
  {
    return("STFEngineNormalize");
  } //  const char const* STFEngineNormalize::name() const

  /*----------------------------------------------------------------------*/

  void STFEngineNormalize::initialize() 
  {
    STFINV_illegal;
  } // void STFEngineNormalize::initialize()

  /*----------------------------------------------------------------------*/

  void STFEngineNormalize::exec() 
  {
    STFINV_illegal;
  } // void STFEngineNormalize::exec()

  /*----------------------------------------------------------------------*/

  void STFEngineNormalize::classhelp(std::ostream& os)
  {
    os << "class STFEngineNormalize (" 
      << STFEngineNormalize::ID << ")\n";
    os << STFEngineNormalize::description << "\n" << std::endl;
    os << "Options and parameters in common for Fourier engines:\n"
      << "NOT YET IMPLEMENTED"
      << std::endl;
    Tbase::classhelp(os);
  } // void STFEngineNormalize::classhelp(std::ostream& os)

} // namespace stfinv

/* ----- END OF stfinvnormalize.cc ----- */
