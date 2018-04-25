#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <signal.h>

#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

#include "libhdr.h"
#include "loki_compress.h"

#define	STDIN	0
#define	STDOUT	1
#define   READ	     0
#define   WRITE	1

static void ignore_handler(int i)
{
  i=i;
  /* Do nothing */
}

int child_open1(const int read_flag,const char *fname,const char *filterprog,const char *arg)
{
  int ppipe[2]={-1,-1},fd= -1,fd1;
  struct stat sbuf;
  struct sigaction s_action;
  int childpid;

  if(read_flag==READ && fname) if(stat(fname,&sbuf)) return fd;
  if(pipe(ppipe)<0)	{
    (void)fprintf(stderr,"child_open(): Can't open pipe\n");
    return fd;
  }
  childpid=fork();
  if(childpid<0)	{
    (void)fprintf(stderr, "child_open(): cannot fork\n");
    return fd;
  }
  if(childpid>0)	{
    if(read_flag==READ)	{
      fd=ppipe[READ];
      if(close(ppipe[WRITE])<0) {
	(void)fprintf(stderr, "child_open(): cannot close pipe\n");
	exit(EXIT_FAILURE);
      }
    } else {
      fd=ppipe[WRITE];
      if(close(ppipe[READ])<0) {
	(void)fprintf(stderr, "child_open(): cannot close pipe\n");
	exit(EXIT_FAILURE);
      }
    }
  } else {
    if(read_flag==READ)	{
      dup2(ppipe[WRITE], STDOUT);
      if(close(ppipe[READ])<0) {
	(void)fprintf(stderr, "child_open(): cannot close pipe\n");
	exit(EXIT_FAILURE);
      }
      if(fname) {
	fd1=open(fname,O_RDONLY,0666);
	if(fd1<0) {
	  (void)fprintf(stderr, "child_open(): cannot open file %s\n",fname);
	  exit(EXIT_FAILURE);
	}
	dup2(fd1,STDIN);
      }
    } else {
      dup2(ppipe[READ], STDIN);
      if(close(ppipe[WRITE])<0) {
	(void)fprintf(stderr, "child_open(): cannot close pipe\n");
	exit(EXIT_FAILURE);
      }
      if(fname) {
	fd1=creat(fname,0666);
	if(fd1<0) {
	  (void)fprintf(stderr, "child_open(): cannot open file %s\n",fname);
	  exit(EXIT_FAILURE);
	}
	dup2(fd1,STDOUT);
      }
    }
    memset(&s_action,0,sizeof(struct sigaction));
    s_action.sa_handler=ignore_handler;
    s_action.sa_flags=0;
    (void)sigaction(SIGHUP,&s_action,0L);
    (void)sigaction(SIGINT,&s_action,0L);
    (void)sigaction(SIGQUIT,&s_action,0L);
    (void)sigaction(SIGPIPE,&s_action,0L);
    if(read_flag==READ) (void)execlp(filterprog,filterprog,arg,(char *)0);
    else (void)execlp(filterprog,filterprog,arg,(char *) 0);
    (void)fprintf(stderr, "child_open(): cannot exec %s\n",filterprog);
    _exit(EXIT_FAILURE);
  }
  return fd;
}

int child_open(const int read_flag,const char *fname,const char *filterprog)
{
  int fd;

  if(read_flag==READ) fd=child_open1(read_flag,fname,filterprog,"-d");
  else fd=child_open1(read_flag,fname,filterprog,0);
  return fd;
}

int child_open_rw(int fd[2],const char *filterprog,char *const argv[]) 
{
  int read_pipe[2]={-1,-1},write_pipe[2]={-1,-1};
  struct sigaction s_action;
  int childpid;

  fd[0]=fd[1]=-1;
  /* Open up a read pipe (from the filter) and a write pipe (to the filter) */
  if(pipe(read_pipe)<0 || pipe(write_pipe)<0)	{
    (void)fprintf(stderr,"child_open_rw(): Can't open pipe\n");
    return -1;
  }
  childpid=fork();
  if(childpid<0)	{
    (void)fprintf(stderr, "child_open_rw(): cannot fork\n");
    return -1;
  }
	
  if(childpid>0)	{
    /* In parent process */

    /* Close write end of read pipe */
    fd[READ]=read_pipe[READ];
    if(close(read_pipe[WRITE])<0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    /* Close read end of write pipe */
    fd[WRITE]=write_pipe[WRITE];
    if(close(write_pipe[READ])<0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
  } else {
    /* In child process */

    /* Duplicate STDOUT to write end of read pipe, and close read end */
    dup2(read_pipe[WRITE],STDOUT);
    if(close(read_pipe[READ])<0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    /* Duplicate STDIN to read end of write pipe, and close write end */
    dup2(write_pipe[READ],STDIN);
    if(close(write_pipe[WRITE])<0) {
      (void)fprintf(stderr, "child_open_rw(): cannot close pipe\n");
      exit(EXIT_FAILURE);
    }
    s_action.sa_handler=ignore_handler;
    s_action.sa_flags=0;
    (void)sigaction(SIGHUP,&s_action,0L);
    (void)sigaction(SIGINT,&s_action,0L);
    (void)sigaction(SIGQUIT,&s_action,0L);
    (void)sigaction(SIGPIPE,&s_action,0L);
    (void)execv(filterprog,argv);
    (void)fprintf(stderr, "child_open_rw(): cannot exec %s\n",filterprog);
    _exit(EXIT_FAILURE);
  }
  return 0;
}

FILE *_open_readfile(const char *fname,int *flag,const struct lk_compress *compress,int chk_flag)
{
  int i,guess=COMPRESS_NONE;
  FILE *fptr;
  unsigned char buf[4];
  char *filter;
  char *prog[]={"gzip","bzip2","zip","compress"};

  errno=0;
  *flag=0;
  if(!fname) return stdin;
  if(!(fptr=fopen(fname,"r"))) {
    if(chk_flag) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for reading (%s)\n",__func__,fname,strerror(errno));
    else {
      fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for reading (%s)\n",__FILE__,__LINE__,__func__,fname,strerror(errno));
      return 0;
    }
  }
  i=(int)fread(buf,(size_t)1,(size_t)4,fptr);
  if(i==4) {
    if(buf[0]==0x1f) {
      if(buf[1]==0x9d) guess=COMPRESS_COMPRESS; /* compress */
      else {
	if(buf[1]==0x8b && buf[2]==0x08) guess=COMPRESS_GZIP; /* gzip */
      }
    } else {
      if(buf[0]=='B' && buf[1]=='Z' && buf[2]=='h' && buf[3]>='0' && buf[3]<='9') guess=COMPRESS_BZIP2; /* bzip2 */
      else {
	if(buf[0]=='P' && buf[1]=='K' && buf[2]==0x03 && buf[3]==0x04) guess=COMPRESS_ZIP; /* zip */
      }
    }
  }
  fclose(fptr);
  if(guess<COMPRESS_NONE) {
    message(DEBUG_MSG,"File appears to have been compressed using %s\n",prog[guess]);
    filter=compress->comp_path[guess][guess==COMPRESS_ZIP?1:0];
    if(filter) {
      *flag=1;
      if(guess==COMPRESS_ZIP) i=child_open1(READ,fname,filter,0);
      else i=child_open1(READ,fname,filter,"-d");
      if(!(fptr=fdopen(i,"r"))) ABT_FUNC("Couldn't fdopen() stream");
      if(errno && errno!=ESPIPE) ABT_FUNC("Unknown IO error\n");
      errno=0;
    } else {
      if(chk_flag) abt(__FILE__,__LINE__,"%s(): File '%s' appears to have been compressed using %s, which is not in the current $PATH\n",__func__,fname,prog[guess]);
      else {
	fprintf(stderr,"[%s:%d] %s(): File '%s' appears to have been compressed using %s, which is not in the current $PATH\n",__FILE__,__LINE__,__func__,fname,strerror(errno));
	fptr=0;
      }
    }
  } else {
    message(DEBUG_MSG,"File is uncompressed\n");
    if(!(fptr=fopen(fname,"r"))) {
      if(chk_flag) abt(__FILE__,__LINE__,"%s(): File Error.  Couldn't open '%s' for reading\n",__func__,fname);
      else fprintf(stderr,"[%s:%d] %s(): File Error.  Couldn't open '%s' for reading (%s)\n",__FILE__,__LINE__,__func__,fname,strerror(errno));
    }
  }
  return fptr;
}
