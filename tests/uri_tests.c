//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.2.0: ExaScale IO library for turbulence simulation restart files
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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

#include <esio/uri.h>

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
    FCT_SUITE_BGN(scheme_prefix_len)
    {
        FCT_TEST_BGN(degenerate)
        {
            fct_chk_eq_int(0, scheme_prefix_len(NULL));
            fct_chk_eq_int(0, scheme_prefix_len(""));
            fct_chk_eq_int(0, scheme_prefix_len("foo"));
            fct_chk_eq_int(0, scheme_prefix_len("foo/bar"));
            fct_chk_eq_int(0, scheme_prefix_len("/foo/bar/baz"));
            fct_chk_eq_int(0, scheme_prefix_len("/foo/bar/baz qux quux"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(empty_path)
        {
            fct_chk_eq_int(2, scheme_prefix_len("s:"));
            fct_chk_eq_int(3, scheme_prefix_len("rs:"));
            fct_chk_eq_int(4, scheme_prefix_len("rst:"));
            fct_chk_eq_int(8, scheme_prefix_len("rst.lne:"));
            fct_chk_eq_int(8, scheme_prefix_len("rst+lne:"));
            fct_chk_eq_int(8, scheme_prefix_len("rst-lne:"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(simple_path)
        {
            fct_chk_eq_int(2, scheme_prefix_len("s:foo"));
            fct_chk_eq_int(3, scheme_prefix_len("rs:bar"));
            fct_chk_eq_int(4, scheme_prefix_len("rst:baz"));
            fct_chk_eq_int(8, scheme_prefix_len("rst.lne:qux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst+lne:quux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst-lne:quuux"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(compound_path)
        {
            fct_chk_eq_int(2, scheme_prefix_len("s:foo/bar"));
            fct_chk_eq_int(3, scheme_prefix_len("rs:bar/baz"));
            fct_chk_eq_int(4, scheme_prefix_len("rst:baz/qux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst.lne:qux/quux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst+lne:quux/quuux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst-lne:quuux/quuuux"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(absolute_path)
        {
            fct_chk_eq_int(2, scheme_prefix_len("s:/foo/bar"));
            fct_chk_eq_int(3, scheme_prefix_len("rs:/bar/baz"));
            fct_chk_eq_int(4, scheme_prefix_len("rst:/baz/qux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst.lne:/qux/quux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst+lne:/quux/quuux"));
            fct_chk_eq_int(8, scheme_prefix_len("rst-lne:/quuux/quuuux"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(trickier)
        {
            // Empty
            fct_chk_eq_int(0, scheme_prefix_len(":/foo/bar"));

            // Does not begin with isalpha
            fct_chk_eq_int(0, scheme_prefix_len("9:/foo/bar"));
            fct_chk_eq_int(0, scheme_prefix_len("+:/foo/bar"));
            fct_chk_eq_int(0, scheme_prefix_len("-:/foo/bar"));
            fct_chk_eq_int(0, scheme_prefix_len(".:/foo/bar"));

            // Colon embedded later in path
            fct_chk_eq_int(0, scheme_prefix_len("/foo:bar/baz"));
            fct_chk_eq_int(0, scheme_prefix_len("/foo/bar:baz qux quux"));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(use_case)
        {
            const char *str1 = "test:ing123";
            fct_chk_eq_str(str1 + scheme_prefix_len(str1), "ing123");

            const char *str2 = "testing123";
            fct_chk_eq_str(str2 + scheme_prefix_len(str2), str2);
        }
        FCT_TEST_END();
    }
    FCT_SUITE_END();
}
FCT_END()
