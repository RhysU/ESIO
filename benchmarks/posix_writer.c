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
#include<string.h>
#include<stdlib.h>

int xfersize=256*1024;         	/* low-level write size */
int blocked=1;			/* Blocked I/O mode? */
int blocksize=4*1024*1024;	/* block write size (bytes) */

/*
  starttime = MPI_Wtime();
    ....  stuff to be timed  ...
    endtime   = MPI_Wtime();
printf("That took %f seconds\n",endtime-starttime)
*/

double starttime;
double endtime;

double writetime;
double donewritetime;

int  fd_global;
FILE *fp_global;
int  lustre_fs=1;  // this is NOT a lustre file system!
int  lustre_stripe_index=-1;    // -1 = default Lustre value
int  lustre_stripe_count=4; 	//  0 = default Lustre value

//extern void posix_write_doubles_1d(int fd,double *val,size_t count,off_t offset, int local_rank);

int IO_COMM=MPI_COMM_WORLD;

void posix_open(char *filename,int local_rank) 
{
  int fd_flag=0;

  starttime= MPI_Wtime();

  // All tasks open the shared solution file and seek to their portion
  // of the file - note that task 0 is first responsible for creating the
  // file, while the other tasks then open for writing.

  MPI_Barrier(IO_COMM);
  
  if(local_rank == 0)
    {

      fd_flag |= O_CREAT | O_RDWR | O_TRUNC | O_NONBLOCK;
            
      if( (fd_global = open(filename,fd_flag,0600)) == -1)
	{
	  printf("** Error (pwrite): raw open failed for %s (rank %i)\n ",
		 filename,local_rank);
	  MPI_Abort(IO_COMM,1); 
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
  endtime= MPI_Wtime();

}

void posix_close(int fd)
{

  //  fsync(fd);
  fclose(fp_global);
  close(fd_global);

}

//                         filename      soln  num_loc_elements  offset      num_local      
void posix_write_double_1d(int fd,double *val,long long count,long long offset, long long local_rank)
{

  //  float filesize_mb;		/* total filesize */
  float num_writes;
  float elem_remain;
  long index;
  int mycount;
  int elem_write,num_blocks,iseg;

  assert(fd  > 0);
  assert(val != NULL);
  assert(count > 0);
  assert(offset >= 0);
  
  writetime = MPI_Wtime();

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
		  fwrite(&val[index],sizeof(double),mycount,fp_global); // here

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
    donewritetime = MPI_Wtime();
}

int  main(int argc,char *argv[])
{

  int num_procs, num_local;
  char mach_name[MPI_MAX_PROCESSOR_NAME];
  int mach_len;

  long long filesize_mb;		/* total filesize */
  long long num_writes, num_blocks;
  long long offset;
  double *soln;
  char str[10];
  char str2[10];
  char filename[100];

  MPI_Init (&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank (MPI_COMM_WORLD, &num_local);
  MPI_Get_processor_name(mach_name,&mach_len);

  // Command-line parsing -- pass total file size to be written
  filesize_mb = (long long)atof(argv[1]);

  // choose file write location and filename
  strcpy(filename, "/intrepid-fs0/users/nick/scratch/data"); // intrepid
  strcat(filename,"/testfile."); //adding filename
  snprintf(str,sizeof(str),"%d",num_procs); // adding number of processors
  strcat(str,"tasks."); 
  snprintf(str2,sizeof(str),"%llu",filesize_mb); // adding filesize
  strcat(filename,str); //combining these strings for total filepath
  strcat(filename,str2);

  //---------------------------------------------------
  // Create some decomposed solution data - simply 
  // create the same on all processors for convenience
  //---------------------------------------------------

  // this is the total number of elements to be written for ALL processors
  long long num_elements = filesize_mb*1024*1024/sizeof(double); 

  // Determine the local file offset
  long long start_index, end_index;
  long long chunk;
  long long num_local_elements;

  chunk = num_elements/num_procs;
  start_index = num_local*chunk;
  end_index   = start_index + chunk - 1;
  offset      = start_index*sizeof(double);

  if(num_local == num_procs - 1)
      end_index = num_elements - 1;

  num_local_elements = end_index - start_index + 1;
  printf("--> Proc %i: file offset = %llu, local indices = (%llu,%llu), num elements = %llu\n", num_local,offset,start_index,end_index,num_local_elements);
  
  if(num_local == 0)
      printf("--> Local Write Size: %llu MB -- ensure avail. core memory less than this\n", sizeof(double)*num_local_elements/(1024*1024));
  
  // set solution size, given total file size and processors
  soln = calloc(num_local_elements,sizeof(double));

  for(int i=0;i<num_local_elements;i++)
    soln[i] = num_local;

  /* POSIX parallel write tests */
  posix_open(filename,num_local);
  posix_write_double_1d(fd_global,&soln[0],num_local_elements,offset,num_local);
  posix_close(fd_global);

  //---------------------
  // Performance Summary
  //---------------------

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

      write_speed_mb = filesize_mb/((endtime-starttime) + (donewritetime-writetime));
      printf("\n");
      printf("  --> start open file                = %5.1f \n",starttime);
      printf("  --> end open file                  = %5.1f \n",endtime  );
      printf("  --> write start                    = %5.1f \n",donewritetime);
      printf("  --> write finish                   = %5.1f \n",writetime  );
      printf("  --> Total File Size                = %llu MB\n",filesize_mb);
      printf("  --> Local File Write Size          = %llu MB\n",filesize_mb/num_procs);
      printf("  --> Total Write Time               = %5.1f Sec\n",(endtime-starttime+donewritetime-writetime));
      printf("  --> Aggregate Write Speed (POSIX)  = %8.3f (MB/sec), %8.3f (GB/sec)\n",
	     write_speed_mb,write_speed_mb/1024);
    }

  
  MPI_Finalize();
  return 0;
  
}





