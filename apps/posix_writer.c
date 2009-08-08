#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/ioctl.h>
#include<fcntl.h>
#include<math.h>
#include<mpi.h>
#include<errno.h>
#include<assert.h>
#include<lustre/lustre_user.h>
#include<hpct.h>

// karl is just testing

int xfersize=256*1024;         	/* low-level write size */
int blocked=1;			/* Blocked I/O mode? */
int blocksize=4*1024*1024;	/* block write size (bytes) */

int  fd_global;
FILE *fp_global;
int  lustre_fs=1;
int  lustre_stripe_index=-1;    // -1 = default Lustre value
int  lustre_stripe_count=4; 	//  0 = default Lustre value

//extern void posix_write_doubles_1d(int fd,double *val,size_t count,off_t offset, int local_rank);

int IO_COMM=MPI_COMM_WORLD;

void posix_open(char *filename,int local_rank) 
{
  int fd_flag=0;

  hpct_timer_begin("posix_open");

  // All tasks open the shared solution file and seek to their portion
  // of the file - note that task 0 is first responsible for creating the
  // file, while the other tasks then open for writing.

  MPI_Barrier(IO_COMM);
  
  if(local_rank == 0)
    {

      fd_flag |= O_CREAT | O_RDWR | O_TRUNC | O_NONBLOCK;

      if(lustre_fs)
	fd_flag |= O_LOV_DELAY_CREATE;
	  
      if( (fd_global = open(filename,fd_flag,0600)) == -1)
	{
	  printf("** Error (pwrite): raw open failed for %s (rank %i)\n ",
		 filename,local_rank);
	  MPI_Abort(IO_COMM,1); 
	}

      if(lustre_fs)       // add specific Lustre striping
	{
	  struct lov_user_md opts = { 0 };
	  opts.lmm_magic          = LOV_USER_MAGIC;
	  opts.lmm_stripe_size    = 4194304;
	  opts.lmm_stripe_offset  = lustre_stripe_index;
	  opts.lmm_stripe_count   = lustre_stripe_count;

	  ioctl(fd_global,LL_IOC_LOV_SETSTRIPE,&opts);

	  //if(ioctl(fd_global,LL_IOC_LOV_SETSTRIPE,&opts) < 0)
	  //perror("Problems with ioctl call");
	}

      MPI_Barrier(IO_COMM);
    }
  else
    {
      MPI_Barrier(IO_COMM);
      fd_flag |= O_RDWR ;
      if( (fd_global = open(filename,fd_flag,0600)) == -1)
	{
	  printf("** Error (pwrite): raw open failed for %s (rank %i)\n",filename,local_rank);
	  MPI_Abort(IO_COMM,1); 
	}
    }

  //  printf("posix_open_complete\n");

  hpct_timer_end("posix_open");

}

void posix_close(int fd)
{

  hpct_timer_begin("posix_close");
  //  fsync(fd);
  fclose(fp_global);
  close(fd_global);
  hpct_timer_end("posix_close");

}

void posix_write_double_1d(int fd,double *val,size_t count,off_t offset, int local_rank)
{

  //  float filesize_mb;		/* total filesize */
  size_t num_writes;
  size_t elem_remain;
  size_t index;
  int mycount;
  int elem_write,num_blocks,iseg;

  assert(fd  > 0);
  assert(val != NULL);
  assert(count > 0);
  assert(offset >= 0);

  hpct_timer_begin("posix_write_1d");

  // note to self: verify mode should make sure the offset and array
  // length provide non-overlapping segments.

  // Seek to the correct file offset for the local processor and
  // associate a file stream for writing

  if (lseek(fd, offset, SEEK_SET) == -1)
    {
      printf("** Error (pwrite): seek failed for proc %i\n",local_rank);
      MPI_Abort(MPI_COMM_WORLD,1);
    }  

  if( (fp_global = fdopen(fd,"wb")) == NULL)
    {
      printf("** Error (pwrite): file open failed.\n");
      MPI_Abort(MPI_COMM_WORLD,1); 
    }

  fseek(fp_global,offset,SEEK_SET);

  // Write some groovy data

    if(!blocked)
      {
	fwrite (val,sizeof(double),count,fp_global);
      }
    else 
      {

	num_writes = floorf(1.0*count*4./blocksize);
	mycount    = xfersize/4;
	index      = offset/4 + 1;
	index      = 0;
	elem_write = 0;
	num_blocks = num_writes;

	while(num_blocks >= 1)
	  {
	    for(iseg=1;iseg<=blocksize/xfersize;iseg++)
	      {

		fwrite(&val[index],sizeof(double),mycount,fp_global);

		index  += mycount;
		offset += xfersize;
	      }
	  
	    num_blocks = num_blocks - 1;
	  }

	elem_remain = count - ( num_writes*blocksize/4 );
      
	while(elem_remain >0)
	  {
	    if(elem_remain >= mycount)
	      {

		fwrite(&val[index],sizeof(double),mycount,fp_global);

		index  += mycount;
		offset += xfersize;
		elem_remain  = elem_remain - mycount;
	      }
	    else
	      {
		fwrite(&val[index],sizeof(double),elem_remain,fp_global);
		break;
	      }
	  }

      }

  hpct_timer_end("posix_write_1d");

}

int  main(int argc,char *argv[])
{

  int num_procs, num_local;
  char mach_name[MPI_MAX_PROCESSOR_NAME];
  int mach_len;

  float filesize_mb;		/* total filesize */
  int num_writes, num_blocks;
  off_t offset;
  double *soln;
  char *filename = "/work/00161/karl/superfly";

  MPI_Init (&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank (MPI_COMM_WORLD, &num_local);
  MPI_Get_processor_name(mach_name,&mach_len);

  // Command-line parsing

  filesize_mb = (float)atof(argv[1]);

  //---------------------------------------------------
  // Create some decomposed solution data - simply 
  // create the same on all processors for convenience
  //---------------------------------------------------

  size_t num_elements = filesize_mb*1024*1024/sizeof(double);

  soln = calloc(num_elements,sizeof(double));

  for(int i=0;i<num_elements;i++)
    soln[i] = i;

  // Determine the local file offset

  int start_index, end_index;
  int chunk;
  int num_local_elements;

  chunk = num_elements/num_procs;

  start_index = num_local*chunk;
  end_index   = start_index + chunk - 1;
  offset      = start_index*sizeof(double);

  if(num_local == num_procs - 1)
    end_index = num_elements - 1;

  num_local_elements = end_index - start_index + 1;

  printf("--> Proc %i: file offset = %10i, local indices = (%10i,%10i), num elements = %10i\n",
	 num_local,offset,start_index,end_index,num_local_elements);

  hpct_timer_init("I/O Testing");

  /* POSIX parallel write tests */

  posix_open(filename,num_local);
  posix_write_double_1d(fd_global,&soln[start_index],num_local_elements,offset,num_local);
  posix_close(fd_global);

  //---------------------
  // Performance Summary
  //---------------------

  hpct_timer_finalize();

  if(num_local == 0)
    {
      printf("\n** Write Performance Summary **\n");
      printf("  --> Master task on %s\n",mach_name);
      printf("  --> Total number of tasks          = %i\n",num_procs);
      printf("  --> Output filename                = %s\n",filename);

      if(lustre_fs)
	{
	  printf("  --> Lustre stripe index            = %i \n",lustre_stripe_index);
	  printf("  --> Lustre stripe count            = %i \n",lustre_stripe_count);
	}

      float write_speed_mb;
      write_speed_mb = filesize_mb/(hpct_timer_elapsedseconds("posix_open") +
				    hpct_timer_elapsedseconds("posix_write_1d") );
      //				    hpct_timer_elapsedseconds("posix_close") );

      printf("\n");
      printf("  --> Total File Size                = %5.1f MB\n",filesize_mb);
      printf("  --> Aggregate Write Speed (POSIX)  = %8.3f (MB/sec), %8.3f (GB/sec)\n",
	     write_speed_mb,write_speed_mb/1024);

      //      hpct_timer_summarize();

    }

  
  MPI_Finalize();
  return 0;
  
}





