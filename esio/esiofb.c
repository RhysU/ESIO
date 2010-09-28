/* these are just simple routines that take the esiof routines and convert */
/* to esiow routines.  esiofb stands for ESIO F.ortran B.indings */

#include "esio.h"

void esiofb_fopen_(char* file, char* trunc)
{
    esio_fopen(file, trunc);
}

void esiofb_fclose_()
{
    esio_fclose();
}

void esiofb_dopen_(int ny, int nx, int nz, char* datasetname, char* precision)
{
    esio_dopen(ny, nx, nz, datasetname, precision);
}

void esiofb_dclose_()
{
    esio_dclose();
}

void esiofb_write_double_(int ny, int stepover, double* data)
{
    esio_write_double(ny, stepover, data);
}

void esiofb_write_double_field_(int ny, int nx, int nz, int nc,
                                int xst, int xsz, int zst, int zsz,
                                double* data, char *filename,
                                char *dataname, char *overwrite,
                                int fname_len, int dname_len, int oname_len)
{
    esio_write_double_field(ny, nx, nz, nc,
                            xst, xsz, zst, zsz,
                            data, filename, dataname, overwrite);
}
