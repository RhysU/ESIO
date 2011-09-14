//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.7: ExaScale IO library for turbulence simulation restart files
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2010, 2011 The PECOS Development Team
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

#include <ctype.h>

#include "uri.h"

int scheme_prefix_len(const char *s)
{
    // RFC 3986 says
    // scheme = ALPHA *( ALPHA / DIGIT / "+" / "-" / "." )

    if (!s || !isalpha(*s)) return 0;

    const char *t = s;
    while (++t) {
        if (*t == ':') return t - s + 1;
        _Bool valid =    isalpha(*t) // Rough guess of frequency
                      || *t == '.'
                      || *t == '-'
                      || isdigit(*t)
                      || *t == '+';
        if (!valid) return 0;
    }

    return 0;
}
