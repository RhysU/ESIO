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

/**
 * Create an output filename template suitable for use with mkstemp(3).
 * The caller must free the returned pointer using free(3).
 *
 * @param dir Directory to be used in the template.
 * @param filename Filename on which the template is based.
 *
 * @return A newly malloc-ed file template suitable for use with mkstemp(3).
 */
char * create_testfiletemplate(const char *dir, const char * const filename);


/**
 * Given a mkstemp(3)-ready testfiletemplate, create a new temporary filename.
 * The caller must free the returned pointer using free(3).
 *
 * @param testfiletemplate Template to use.
 *
 * @return A new temporary filename based on testfiletemplate.
 *         On failure, returns NULL.
 */
char * create_testfilename(const char * const testfiletemplate);
