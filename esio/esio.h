/* Public C API for ESIO */

#ifndef __ESIO_ESIO_H
#define __ESIO_ESIO_H

#include <mpi.h>

typedef struct esio_state_s *esio_state;

esio_state esio_init(MPI_Comm comm);
void esio_finalize(esio_state s);

int esio_file_create(esio_state s, const char *file, int overwrite);
int esio_file_open(esio_state s, const char *file, int readwrite);
int esio_file_close(esio_state s);

int esio_field_write_double(esio_state s,
                            const char* name,
                            double *data,
                            int na, int ast, int asz,
                            int nb, int bst, int bsz,
                            int nc, int cst, int csz);
int esio_field_write_float(esio_state s,
                           const char* name,
                           float *data,
                           int na, int ast, int asz,
                           int nb, int bst, int bsz,
                           int nc, int cst, int csz);

#endif /* __ESIO_ESIO_H */
