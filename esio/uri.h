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

#ifndef ESIO_URI_H
#define ESIO_URI_H

//****************************************************************
// INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL
//****************************************************************

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Determine the length of the optional URI scheme prefix on \c s, including
 * any trailing colon (<tt>:</tt>).  URI schemes follow the syntax given in <a
 * href="http://www.ietf.org/rfc/rfc3986.txt">RFC 3986</a> section 3 "Syntax
 * components".
 *
 * The length of the prefix <em>includes</em> the trailing colon.  Given string
 * \c s, the following facts hold:
 * <ol>
 * <li><tt>uri_scheme_len(s)</tt> is true iff \c s has a URI prefix.</li>
 * <li><tt>s + uri_scheme_len(s)</tt> points to the non-scheme portion of
 *     the string</li>
 * <li><tt>uri_scheme_len(NULL) == 0</tt>, always.
 * </ol>
 *
 * @param s String which may optionally have a leading URI scheme.
 *
 * @return Length of the URI scheme plus the trailing colon.
 */
int scheme_prefix_len(const char *s);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ESIO_URI_H */
