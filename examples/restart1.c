#include <stdlib.h>
#include <mpi.h>
#include <esio/esio.h>

int main(int argc, char *argv[])
{
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /* Create restart file template */
    esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
    esio_file_create(h, "template.h5", 1 /* overwrite */);
    esio_string_set(h, "/", "program", argv[0]);
    esio_file_close(h);

    /* Establish the parallel decomposition */
    esio_line_establish(h, 2*world_size, 2*world_rank, 2);

    /* Main simulation loop */
    double example[2];
    for (int i = 0; i < 10; ++i) {

        /* Update simulation state */
        example[0] = i * world_rank;
        example[1] = i * world_rank + 1;

        /* Create an uncommitted restart file based on the template */
        esio_file_clone(h, "template.h5", "uncommitted.h5", 1 /*overwrite*/);

        /* Save the current simulation state in the uncommitted restart */
        esio_line_write_double(h, "example", example, 1, NULL);

        /* Commit the restart file.  Up to 5 older restarts are retained. */
        esio_file_close_restart(h, "committed#.h5", 5 /*retain_count*/);
    }

    esio_handle_finalize(h);

    return 0;
}
