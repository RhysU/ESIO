#include <stdlib.h>
#include <mpi.h>
#include <esio/esio.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
    esio_file_create(h, "data.h5", 1 /* overwrite */);

    esio_string_set(h, "program", argv[0]);

    int version = 1;
    esio_attribute_write_int(h, "version", &version);

    double example[2] = { 2.0 * world_rank, 2.0 * world_rank + 1 };
    esio_line_establish(h, 2*world_size, 2*world_rank, 2);
    esio_line_write_double(h, "example", example, 1);

    esio_file_close(h);
    esio_handle_finalize(h);

    return 0;
}
