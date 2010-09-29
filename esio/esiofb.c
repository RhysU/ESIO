/* these are just simple routines that take the esiof routines and convert */
/* to esiow routines.  esiofb stands for ESIO F.ortran B.indings */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "esio.h"

/*
 * See GNU autoconf documentation for AC_FC_WRAPPERS for info on FC_FUNC_:
 * http://www.gnu.org/software/hello/manual/autoconf/Fortran-Compiler.html
 */
#ifndef FC_FUNC_
#error FC_FUNC_ not available from autotools
#endif

#ifdef __INTEL_COMPILER
/* remark #1418: external function definition with no prior declaration */
#pragma warning(disable:1418)
#endif

void FC_FUNC_(esiofb_fopen,ESIOFB_FOPEN)
    (char* file,
     char* trunc)
{
    esio_fopen(file, trunc);
}

void FC_FUNC_(esiofb_fclose,ESIOFB_FCLOSE)
    ()
{
    esio_fclose();
}

void FC_FUNC_(esiofb_dopen,ESIOFB_DOPEN)
    (int *ny,
     int *nx,
     int *nz,
     char *datasetname,
     char *precision)
{
    esio_dopen(*ny, *nx, *nz, datasetname, precision);
}

void FC_FUNC_(esiofb_dclose,ESIOFB_DCLOSE)
    ()
{
    esio_dclose();
}

void FC_FUNC_(esiofb_write_double,ESIOFB_WRITE_DOUBLE)
    (int *ny,
     int *stepover,
     double *data)
{
    esio_write_double(*ny, *stepover, data);
}

void FC_FUNC_(esiofb_write_double_field,ESIOFB_WRITE_DOUBLE_FIELD)
    (int *ny,
     int *nx,
     int *nz,
     int *nc,
     int *xst,
     int *xsz,
     int *zst,
     int *zsz,
     double* data,
     char *filename,
     char *dataname,
     char *overwrite,
     int fname_len,
     int dname_len,
     int oname_len)
{
    (void) fname_len; /* Unused */
    (void) dname_len; /* Unused */
    (void) oname_len; /* Unused */
    esio_write_double_field(*ny, *nx, *nz, *nc,
                            *xst, *xsz, *zst, *zsz,
                            data, filename, dataname, overwrite);
}
