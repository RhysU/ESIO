/* these are the esio writing c routines for hdf5 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <mpi.h>
#include "esio.h"

#define RANK (2)                               /* Dimensionality */
static const MPI_Comm comm =  MPI_COMM_WORLD;  /* Communicator */
static       hid_t dset_id = -1;               /* file identifier */
static       hid_t file_id = -1;               /* dataset identifier */

int esio_fopen(char* file, char* trunc)
{
    hid_t plist_id;                 /* property list identifier */

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

    if (file_id != -1) /* another open file exists-- abort */
    {
        printf("\nESIO FATAL ERROR: "
               "Cannot open new file -- previous file not closed!\n");
        printf("ESIO Does not support multiple open files.\n");
        MPI_Abort(comm, 1);
    }

    /*
     * Create a new file collectively
     * OR
     * Delete old data and update values
     */
    if (strncmp(trunc, "o", 1) == 0) /* starts with n */
    {
        file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        if (file_id < 0)
        {
            printf("\nESIO FATAL ERROR: "
                   "Unable to create file\n");
            MPI_Abort(comm, 1);
        }
    }
    else /* do not overwrite old dataset -- will fail if file already exists */
    {
        file_id = H5Fcreate(file, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
        if (file_id < 0)
        {
            printf("\nESIO FATAL ERROR: "
                   "File already exists\n");
            MPI_Abort(comm, 1);
        }
    }

    H5Pclose(plist_id);   /* release property list identifier */

    return 0;
}

int esio_fclose()
{

    int error = H5Fclose(file_id);

    if (error < 0)
    {
        printf("\nESIO FATAL ERROR: "
               "File unable to close\n");
        MPI_Abort(comm, 1);
    }

    file_id = -1; /* reset unique file identifier */
    return 0;
}

int esio_dopen(int ny, int nx, int nz, char* datasetname, char* precision)
{
    hid_t   filespace;                /* file identifier */
    hsize_t dimsf[RANK];              /* dataset dimensions */

    assert(ny > 0);
    assert(nx > 0);
    assert(nz > 0);

    if (dset_id != -1) /* another open dataset exists -- abort */
    {
        printf("\nESIO FATAL ERROR: "
               "Cannot open new dataset -- previous not closed!\n");
        printf("ESIO Does not support multiple open datasets.\n");
        MPI_Abort(comm, 1);
    }

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = nx * nz;
    dimsf[1] = ny;
    filespace = H5Screate_simple(RANK, dimsf, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    if (strncmp(precision, "d", 1) == 0)
    {
        /* write doubles */
        dset_id = H5Dcreate1(file_id, datasetname, H5T_NATIVE_DOUBLE,
                             filespace, H5P_DEFAULT);
    }
    else if ((strncmp(precision, "f", 1) || strncmp(precision, "s", 1)) == 0)
    {
        /* write single precision */
        dset_id = H5Dcreate1(file_id, datasetname, H5T_NATIVE_FLOAT,
                             filespace, H5P_DEFAULT);
    }
    else
    {
        printf("\nESIO FATAL ERROR: "
                "Cannot open dataset -- precision not properly specified\n");
        printf("User Must specify either double or single precision\n");
        MPI_Abort(comm, 1);
    }

    H5Sclose(filespace);

    return 0;

}

int esio_dclose()
{

    int error = H5Dclose(dset_id);
    if (error < 0)
    {
        printf("\nESIO FATAL ERROR: "
               "File unable to close\n");
        MPI_Abort(comm, 1);
    }
    dset_id = -1; /* reset */
    return 0;
}



/*
 * Each process defines dataset in memory and writes it to the hyperslab
 * in the file.
 */

static int esio_write_double(int ny, int zxstepover, double* data)
{

    hid_t    filespace;      /* file and memory dataspace identifiers */
    hid_t    memspace;       /* memory dataspace identifier */
    hid_t    plist_id;       /* property list identifier */
    herr_t  status;

    /* hyperslab selection parameters */
     hsize_t count[RANK];
     hsize_t offset[RANK];
     hsize_t stride[RANK];
     count[0]  = 1;          /* # of values to write */
     count[1]  = ny;
     offset[0] = zxstepover; /* offset               */
     offset[1] = 0;
     stride[0] = 1;
     stride[1] = 1;          /* dont skip values */

    /*
     * Initialize data buffer
     */
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                        offset, stride, count, NULL);

    /*
     * Create property list for collective dataset write
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,
                      filespace, plist_id, data);
    if (status < 0)
    {
        /* user wants to know if data isn't being written properly */
        printf("\nESIO FATAL ERROR: Write Failed\n");
        MPI_Abort(comm, 1);
    }

    /*
     * Close/release resources
     */
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    return 0;

}

int esio_write_double_field(int ny, int nx, int nz, int nc,
                            int xst, int xsz,
                            int zst, int zsz, double* data,
                            char *filename, char *dataname, char *overwrite)
{
#ifdef TIMERS
    double start, end, timer;
#endif
    int mpi_rank;
    int i, j;

#ifdef TIMERS
    start = MPI_Wtime();                                  /* start timers */
#endif

    MPI_Comm_rank(comm, &mpi_rank);
    esio_fopen(filename, overwrite);                      /* create file */
    esio_dopen(ny, nx, nz, dataname, "double");           /* create dataset */

    --zst; /* convert to zero-indexed value */            /* FIXME Unify */
    --xst; /* convert to zero-indexed value */            /* FIXME Unify */
    for (i = 0; i < xsz; i++)
        for (j = 0; j < zsz; j++)
        {
            const int zxstepover = (j + zst) + (i + xst) * nz;
            const int loffset    = j*ny + i*ny*zsz;
            esio_write_double(ny, zxstepover, data + loffset);
        }

    esio_dclose();                                        /* close dataset */
    esio_fclose();                                        /* close file */

#ifdef TIMERS
    end = MPI_Wtime();                                    /* end timers */
#endif

    /* finalize and report */
#ifdef TIMERS
    if (mpi_rank == 0)
    {
        timer = end - start;
        printf("\nfinalizing\n");
        printf("total time was: %g\n", timer);
        printf("wrote: %g MB\n",
               (nc*sizeof(double)*ny*nz*nx) / (1024.*1024.));
        printf("write speed was: %g GB/s\n",
               ((nc*sizeof(double)*ny*nz*nx) / (1024.*1024.*1024.)) / (timer));
/*         printf("wrote: %g MB\n", */
/*                (nc*sizeof(double)*ny*nz*nx/2.) / (1024.*1024.)); */
    }
#endif

    return 0;
}
