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

#ifndef __ESIO_RESTART_RENAME_H
#define __ESIO_RESTART_RENAME_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Increment any index number found in \c name when it matches \c tmpl.
 * A match requires that the strings are character-by-character identical
 * up to the first \c # (hash sign) and after the last \c # (hash sign)
 * in \c tmpl.  The string \c name must contain one or more decimal
 * digits within the \c # region within \c tmpl.  The current index number
 * is defined to be the decimal value parsed from the \c # region within
 * \c name.  The incremented index number is one plus the current index.
 *
 * @param tmpl   Template to use.
 * @param name   Name to use.
 * @param errval Value to return on usage error due to a bad template
 *               or an index overflow condition.
 *
 * @return The incremented index number on success or zero if \c name
 *         does not match \c tmpl.  Note that a successful, matching result
 *         will be strictly positive whenever \c errval < 0.
 */
int restart_nextindex(const char *tmpl,
                      const char *name,
                      const int errval);

int restart_rename(const char *src_filename,
                   const char *dst_template,
                   int keep_howmany);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __ESIO_RESTART_RENAME_H */
