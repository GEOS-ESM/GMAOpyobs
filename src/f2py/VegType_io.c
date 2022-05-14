/*

Simple C based I/O to be called from Fortran. Loosely based on similar
routines from RAMS.

 */

#include <stdio.h>
#include <stdlib.h>

static FILE *file;

int vegopen(char *filename,char *mode)
{
  file=fopen(filename,mode);
  if ( file==NULL ) {
    printf ("vegopen: cannot open file <%s>\n", filename);
    exit(1);
  }
  return(0);
}
int vegclose()
{
  return (int) fclose(file);
}
int vegchar_(int *pos, int *size, int *a)
{
  int rc;
  // printf("Now we are in vegchar, pos, size = %d, %d\n", *pos, *size);
   rc=fseek(file,*pos,0);
   fread(a,1,*size,file);
  return(rc);
}

/* Fortran prototypes */
int vegopen_(char *filename,char *mode) { return vegopen(filename,mode); }
int vegclose_()                         { return vegclose(); }
// int vegchar_(int *pos,int *size,int *a) {return vegchar(pos,size,a);}

int VEGOPEN(char *filename,char *mode) { return vegopen(filename,mode); }
int VEGCLOSE()                         { return vegclose(); }
// int VEGCHAR(int *pos,int *size,int *a) {return vegchar(pos,size,a);}

int VEGOPEN_(char *filename,char *mode) { return vegopen(filename,mode); }
int VEGCLOSE_()                         { return vegclose(); }
// int VEGCHAR_(int *pos,int *size,int *a) {return vegchar(pos,size,a);}






