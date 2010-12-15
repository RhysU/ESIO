//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.0.1: ExaScale IO library for turbulence simulation restart files
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <esio/restart-rename.h>

#include "testutils.h"

// Include FCTX and silence useless warnings
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:981)
#endif
#include "fct.h"
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

FCT_BGN()
{
    FCT_FIXTURE_SUITE_BGN(restart_helpers)
    {
        FCT_SETUP_BGN()
        {
            // NOP
        }
        FCT_SETUP_END();

        FCT_TEARDOWN_BGN()
        {
            // NOP
        }
        FCT_TEARDOWN_END();

        FCT_TEST_BGN(nextindex_successful)
        {
            int (* const impl)(const char *, const char *, const int)
                = &restart_nextindex; // Shorthand

            // Middle, one hash
            fct_chk_eq_int(  1, impl("abc#def", "abc0def"  , -1));
            fct_chk_eq_int(  1, impl("abc#def", "abc00def" , -1));
            fct_chk_eq_int( 10, impl("abc#def", "abc9def"  , -1));
            fct_chk_eq_int( 10, impl("abc#def", "abc09def" , -1));
            fct_chk_eq_int( 10, impl("abc#def", "abc009def", -1));
            fct_chk_eq_int(124, impl("abc#def", "abc123def", -1));

            // Leading, one hash
            fct_chk_eq_int(  1, impl("#def", "0def"  , -1));
            fct_chk_eq_int(  1, impl("#def", "00def" , -1));
            fct_chk_eq_int( 10, impl("#def", "9def"  , -1));
            fct_chk_eq_int( 10, impl("#def", "09def" , -1));
            fct_chk_eq_int( 10, impl("#def", "009def", -1));
            fct_chk_eq_int(124, impl("#def", "123def", -1));

            // Trailing, one hash
            fct_chk_eq_int(  1, impl("abc#", "abc0"  , -1));
            fct_chk_eq_int(  1, impl("abc#", "abc00" , -1));
            fct_chk_eq_int( 10, impl("abc#", "abc9"  , -1));
            fct_chk_eq_int( 10, impl("abc#", "abc09" , -1));
            fct_chk_eq_int( 10, impl("abc#", "abc009", -1));
            fct_chk_eq_int(124, impl("abc#", "abc123", -1));

            // Middle, multiple hashes
            fct_chk_eq_int(  1, impl("abc###def", "abc0def"  , -1));
            fct_chk_eq_int(  1, impl("abc###def", "abc00def" , -1));
            fct_chk_eq_int( 10, impl("abc###def", "abc9def"  , -1));
            fct_chk_eq_int( 10, impl("abc###def", "abc09def" , -1));
            fct_chk_eq_int( 10, impl("abc###def", "abc009def", -1));
            fct_chk_eq_int(124, impl("abc###def", "abc123def", -1));

            // Leading, multiple hashes
            fct_chk_eq_int(  1, impl("###def", "0def"  , -1));
            fct_chk_eq_int(  1, impl("###def", "00def" , -1));
            fct_chk_eq_int( 10, impl("###def", "9def"  , -1));
            fct_chk_eq_int( 10, impl("###def", "09def" , -1));
            fct_chk_eq_int( 10, impl("###def", "009def", -1));
            fct_chk_eq_int(124, impl("###def", "123def", -1));

            // Trailing, multiple hashes
            fct_chk_eq_int(  1, impl("abc###", "abc0"  , -1));
            fct_chk_eq_int(  1, impl("abc###", "abc00" , -1));
            fct_chk_eq_int( 10, impl("abc###", "abc9"  , -1));
            fct_chk_eq_int( 10, impl("abc###", "abc09" , -1));
            fct_chk_eq_int( 10, impl("abc###", "abc009", -1));
            fct_chk_eq_int(124, impl("abc###", "abc123", -1));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(nextindex_nomatch)
        {
            int (* const impl)(const char *, const char *, const int)
                = &restart_nextindex; // Shorthand

            // Empty name, one hash
            fct_chk_eq_int(  0, impl("abc#def", "", -1));

            // Leading difference, one hash
            fct_chk_eq_int(  0, impl("abc#def", "ABC0def", -1));
            fct_chk_eq_int(  0, impl("abc#def", "ab0def" , -1));
            fct_chk_eq_int(  0, impl("abc#def", "bc0def" , -1));
            fct_chk_eq_int(  0, impl("abc#def", "0def"   , -1));

            // Trailing difference, one hash
            fct_chk_eq_int(  0, impl("abc#def", "abc0DEF", -1));
            fct_chk_eq_int(  0, impl("abc#def", "abc0de" , -1));
            fct_chk_eq_int(  0, impl("abc#def", "abc0ef" , -1));
            fct_chk_eq_int(  0, impl("abc#def", "abc0"   , -1));

            // No current index and/or leading signs, one hash
            fct_chk_eq_int(  0, impl("abc#def", "abcdef",   -1));
            fct_chk_eq_int(  0, impl("abc#def", "abc+1def", -1));
            fct_chk_eq_int(  0, impl("abc#def", "abc-1def", -1));

            // Empty name, multiple hashes
            fct_chk_eq_int(  0, impl("abc###def", "", -1));

            // Leading difference, multiple hashes
            fct_chk_eq_int(  0, impl("abc###def", "ABC0def", -1));
            fct_chk_eq_int(  0, impl("abc###def", "ab0def" , -1));
            fct_chk_eq_int(  0, impl("abc###def", "bc0def" , -1));
            fct_chk_eq_int(  0, impl("abc###def", "0def"   , -1));

            // Trailing difference, multiple hashes
            fct_chk_eq_int(  0, impl("abc###def", "abc0DEF", -1));
            fct_chk_eq_int(  0, impl("abc###def", "abc0de" , -1));
            fct_chk_eq_int(  0, impl("abc###def", "abc0ef" , -1));
            fct_chk_eq_int(  0, impl("abc###def", "abc0"   , -1));

            // No current index and/or leading signs, multiple hashes
            fct_chk_eq_int(  0, impl("abc###def", "abcdef",   -1));
            fct_chk_eq_int(  0, impl("abc###def", "abc+1def", -1));
            fct_chk_eq_int(  0, impl("abc###def", "abc-1def", -1));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(nextindex_usageerror)
        {
            int (* const impl)(const char *, const char *, const int)
                = &restart_nextindex; // Shorthand

            // Invalid template with default error value
            fct_chk_eq_int(-1, impl("abcdef",      "abcdef",  -1));
            fct_chk_eq_int(-1, impl("abc#!#def",   "abc0def", -1));
            fct_chk_eq_int(-1, impl("abc##!#def",  "abc0def", -1));
            fct_chk_eq_int(-1, impl("abc#!##def",  "abc0def", -1));
            fct_chk_eq_int(-1, impl("abc##!##def", "abc0def", -1));

            // Invalid template with weird error value
            fct_chk_eq_int(567, impl("abcdef",      "abcdef",  567));
            fct_chk_eq_int(567, impl("abc#!#def",   "abc0def", 567));
            fct_chk_eq_int(567, impl("abc##!#def",  "abc0def", 567));
            fct_chk_eq_int(567, impl("abc#!##def",  "abc0def", 567));
            fct_chk_eq_int(567, impl("abc##!##def", "abc0def", 567));

            // Invalid template with zero error value
            fct_chk_eq_int(0, impl("abcdef",      "abcdef",  0));
            fct_chk_eq_int(0, impl("abc#!#def",   "abc0def", 0));
            fct_chk_eq_int(0, impl("abc##!#def",  "abc0def", 0));
            fct_chk_eq_int(0, impl("abc#!##def",  "abc0def", 0));
            fct_chk_eq_int(0, impl("abc##!##def", "abc0def", 0));
        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
