program main

    use mpi
    use esio

    implicit none

    integer             :: ierr, world_size, world_rank, version
    type(esio_handle)   :: h
    character(len = 255):: str
    double precision    :: example(2)

    call MPI_Init (ierr)
    call MPI_Comm_size (MPI_COMM_WORLD, world_size, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, world_rank, ierr)

    call esio_handle_initialize (h, MPI_COMM_WORLD)
    call esio_file_create (h, "data.h5", .true.)

    call get_command_argument (0, str)
    call esio_string_set (h, "program", trim(str))

    version = 1
    call esio_attribute_write_integer (h, "version", version)

    example = [2d0 * world_rank, 2d0 * world_rank + 1]
    call esio_line_establish (h, 2*world_size, 2*world_rank + 1, 2)
    call esio_line_write_double (h, "example", example, 1)

    call esio_file_close (h)
    call esio_handle_finalize (h)

    call MPI_Finalize (ierr)

end program main
