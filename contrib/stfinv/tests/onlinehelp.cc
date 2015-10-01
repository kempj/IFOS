/*! \file onlinehelp.cc
 * \brief print online help
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: onlinehelp.cc 3924 2011-05-09 14:52:03Z tforb $
 * \author Thomas Forbriger
 * \date 09/05/2011
 * 
 * print online help
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
 *  - 09/05/2011   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define ONLINEHELP_VERSION \
  "ONLINEHELP   V1.0   print online help"
#define ONLINEHELP_CVSID \
  "$Id: onlinehelp.cc 3924 2011-05-09 14:52:03Z tforb $"

#include <iostream>
#include <tfxx/commandline.h>
#include <stfinv/stfinvany.h>

using std::cout;
using std::cerr;
using std::endl;

int main(int iargc, char* argv[])
{

  // define usage information
  char usage_text[]=
  {
    ONLINEHELP_VERSION "\n"
    "usage: onlinehelp" "\n"
    "   or: onlinehelp --help|-h" "\n"
  };

  // define full help text
  char help_text[]=
  {
    ONLINEHELP_CVSID
  };

  // define commandline options
  using namespace tfxx::cmdline;
  static Declare options[]= 
  {
    // 0: print help
    {"help",arg_no,"-"},
    // 1: verbose mode
    {"v",arg_no,"-"},
    {NULL}
  };

  // no arguments? print usage...
  if (iargc<2) 
  {
    cerr << usage_text << endl;
    stfinv::STFEngine::help(cerr);
    exit(0);
  }

  // collect options from commandline
  Commandline cmdline(iargc, argv, options);

  // help requested? print full help text...
  if (cmdline.optset(0))
  {
    cerr << usage_text << endl;
    cerr << help_text << endl;
    stfinv::STFEngine::help(cerr);
    stfinv::help(cerr);
    exit(0);
  }

  // dummy operation: print option settings
  for (int iopt=0; iopt<2; iopt++)
  {
    cout << "option: '" << options[iopt].opt_string << "'" << endl;
    if (cmdline.optset(iopt)) {  cout << "  option was set"; }
    else { cout << "option was not set"; }
    cout << endl;
    cout << "  argument (string): '" << cmdline.string_arg(iopt) << "'" << endl;
    cout << "     argument (int): '" << cmdline.int_arg(iopt) << "'" << endl;
    cout << "    argument (long): '" << cmdline.long_arg(iopt) << "'" << endl;
    cout << "   argument (float): '" << cmdline.float_arg(iopt) << "'" << endl;
    cout << "  argument (double): '" << cmdline.double_arg(iopt) << "'" << endl;
    cout << "    argument (bool): '";
    if (cmdline.bool_arg(iopt))
    { cout << "true"; } else { cout << "false"; }
    cout << "'" << endl;
  }
  while (cmdline.extra()) { cout << cmdline.next() << endl; }

  // dummy operation: print rest of command line
  while (cmdline.extra()) { cout << cmdline.next() << endl; }
}

/* ----- END OF onlinehelp.cc ----- */
