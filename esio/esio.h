/* Public C API for ESIO */

#ifndef __ESIO_ESIO_H
#define __ESIO_ESIO_H

int esio_fopen(char* file, char* trunc);
int esio_fclose();
int esio_dopen(int ny, int nx, int nz, char* datasetname, char* precision);
int esio_dclose();
int esio_write_double(int ny, int stepover, double* data);
int esio_write_double_field(int ny, int nx, int nz, int nc,
                            int xst, int xsz,
                            int zst, int zsz, double* data,
                            char *filename, char *dataname, char *overwrite);

#endif /* __ESIO_ESIO_H */
