#include <stdlib.h>
#include <mpi.h>
#include <esio/esio.h>

int main(int argc, char *argv[])
{
    // Initialize MPI and arrange for its finalization
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Create restart file template filled with metadata on a subset of ranks
    if (world_rank == 0) {
        esio_handle s = esio_handle_initialize(MPI_COMM_SELF);
        esio_file_create(s, "template.h5", 1 /* overwrite */);
        esio_string_set(s, "/", "program", argv[0]);
        esio_string_set(s, "/", "built",   __DATE__ " " __TIME__);
        esio_handle_finalize(s);
    }
    MPI_Barrier(MPI_COMM_WORLD);  // Prevent template create/clone race

    // Prepare a handle informed about the parallel domain decomposition
    esio_handle w = esio_handle_initialize(MPI_COMM_WORLD);
    esio_line_establish(w, 2*world_size, 2*world_rank, 2);

    // Main simulation loop would customarily be driven by time advancement
    double state[2];
    for (int i = 0; i < 10; ++i) {

        // Pretend to update simulation state for demonstration purposes
        state[0] = i * world_rank;
        state[1] = i * world_rank + 1;

        // Create an uncommitted restart file based on the template
        esio_file_clone(w, "template.h5", "uncommitted.h5", 1 /*overwrite*/);

        // Save the current simulation state in the uncommitted restart
        esio_line_write_double(w, "state", state, 1, "A comment");

        // Commit the restart file.  Here, up to 5 older restarts are retained
        esio_file_close_restart(w, "committed#.h5", 5 /*retain_count*/);
    }

    esio_handle_finalize(w);
}
