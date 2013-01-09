//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3.0 of the License, or
// (at your option) any later version.
//
// ESIO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
//
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "uri.h"

#include <ctype.h>

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
