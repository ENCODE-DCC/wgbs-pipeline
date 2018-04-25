#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <sys/un.h>
#if HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <errno.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/stat.h>
#if HAVE_SYS_SYSTEMINFO_H
#include <sys/systeminfo.h>
#endif

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 64
#endif

#include "utils.h"
#include "loki.h"
#include "count_dbr.h"
#include "snprintf.h"
#include "loki_monitor.h"

#define MAX_CONN 64
#define MAX_CHILD 8

int child_alive,lmon_shm_id=-1;

static pid_t monitor_pid,child_id[MAX_CHILD];
static char *monitor_sock_addr=".loki_monitor";
static char *loki_pid_file=".loki_pid";
static int n_child;

void reaper(int i)
{
	pid_t id;
	
	do {
		id=waitpid(-1,&i,WNOHANG|WUNTRACED);
		if(id==monitor_pid) {
			if(WIFSTOPPED(i)) {
				(void)fprintf(stderr,"Child (%d) has stopped - sending SIGCONT\n",(int)monitor_pid);
				(void)kill(id,SIGCONT);
			} else {
				(void)unlink(monitor_sock_addr);
				(void)fprintf(stderr,"Child (%d) has died: ",(int)monitor_pid);
				if(WIFEXITED(i)) (void)fprintf(stderr,"exited with status %d\n",WEXITSTATUS(i));
				else (void)fprintf(stderr,"Terminated by signal %d\n",WTERMSIG(i));
				child_alive=0;
			}
		}
	} while(id>0);
}

void reaper1(int i)
{
	pid_t id,j;
	
	do {
		id=waitpid(-1,&i,WNOHANG|WUNTRACED);
		if(id>0 && WIFEXITED(i)) {
			for(j=0;j<n_child;j++) if(child_id[j]==id) {
				child_id[j]=child_id[--n_child];
				break;
			}
		}
	} while(id>0);
}

static int get_num(char *p,int *v,int n,int *fg)
{
	int i;
	char *p1;
	
	for(i=0;i<n;i++) {
		if(i && *(p++)!=',') break;
		v[i]=(int)strtol(p,&p1,10);
		fg[i]=(p==p1)?0:1;
		p=p1;
	}
	return i;
}

static void monitor(int s)
{
	int i,j,k,nfd,bp,bp1,sock,socklist[MAX_CONN],nsock=0,buflen=0,er;
	int v[5],mask[5];
	pid_t par_id,child;
	char buf[256],buf1[8];
	socklen_t len;
	fd_set readfds,tmpfds;
	struct timeval timeout;
	double z,z1,z2,z3,utime[LMON_WIN_SIZE+1],tps;
	time_t tm[LMON_WIN_SIZE+1];
	int iter[LMON_WIN_SIZE+1];
	struct lmon_param *lpar;
	struct sockaddr_un addr;
	struct sigaction s_action;
	char *cmdlist[]={"LMONSTAT","LMONINFO","LSETINFO","LMON_PID","DBR_INFO",
		  "LSIGSTOP","LSIGCONT","LSIGQUIT","LSIGKILL",0};
	
	n_child=0;
	s_action.sa_handler=ignore_handler;
	s_action.sa_flags=0;
	(void)sigaction(SIGPIPE,&s_action,0L);
	s_action.sa_handler=reaper1;
	s_action.sa_flags=0;
	(void)sigaction(SIGCHLD,&s_action,0L);
	lpar=shmat(lmon_shm_id,0,0);
	if(lpar==(void *)-1) {
		par_id=getppid();
		if(par_id!=1) {
			(void)kill(par_id,SIGQUIT);
			perror("Failed to attach shared memory segment");
			ABT_FUNC(AbMsg);
		}
		(void)unlink(monitor_sock_addr);
		(void)unlink(loki_pid_file);
		if(lpar && lpar->dbr_shm_id>=0) free_dbr_count();
		(void)shmctl(lmon_shm_id,IPC_RMID,0);
		_exit(EXIT_SUCCESS);
	}
	bp=bp1=0;
	tps=(double)sysconf(_SC_CLK_TCK);
	for(;;) {
		FD_ZERO(&readfds);
		FD_SET(s,&readfds);
		for(i=0;i<nsock;i++) FD_SET(socklist[i],&readfds);
		/* Set timeout to 1 second to check 
		 * if parent is still alive */
		timeout.tv_sec=1;
		timeout.tv_usec=0;
		nfd=select(FD_SETSIZE,&readfds,0,0,&timeout);
		/* If there is an error, or if the parent's pid is 1 
		 * (which means the original parent has died), then 
		 * cleanup and exit */
		par_id=getppid();
		if((nfd<0 && errno!=EINTR) || par_id==1) {
			for(i=0;i<nsock;i++) {
				(void)shutdown(socklist[i],SHUT_RDWR);
				(void)close(socklist[i]);
			}
			(void)unlink(monitor_sock_addr);
			if(par_id==1) (void)unlink(loki_pid_file);
			if(lpar && lpar->dbr_shm_id>=0) free_dbr_count();
			(void)shmctl(lmon_shm_id,IPC_RMID,0);
			if(nfd<0 && errno!=EINTR) {
				if(par_id!=1) (void)kill(par_id,SIGQUIT);
				perror("Select failed");
				ABT_FUNC(AbMsg);
			}
			_exit(EXIT_SUCCESS);
		}
		if(nfd<0 || lpar->magic!=LMON_MAGIC) continue;
		if(!nfd) continue;		
		iter[bp]=lpar->it;
		utime[bp]=(double)lpar->utime/tps;
		tm[bp]=time(0);
		j=iter[bp]-iter[bp1];
		if(j) {
			z=difftime(tm[bp],tm[bp1])/(double)j;
			z1=(utime[bp]-utime[bp1])/(double)j;
		} else z=z1=0.0;
		z2=lpar->extra_time+difftime(tm[bp],start_time);
		z3=lpar->extra_utime+utime[bp];
		bp++;
		if(bp>LMON_WIN_SIZE) bp=0;
		if(bp1>=bp) {
			bp1++;
			if(bp1>LMON_WIN_SIZE) bp1=0;
		}
		timeout.tv_sec=0;
		timeout.tv_usec=0;
		if(FD_ISSET(s,&readfds)) {
			len=sizeof(addr);
			if((sock=accept(s,(struct sockaddr *)&addr,&len))<0) {
				(void)fprintf(stderr,"[%s:%d] %s():",__FILE__,__LINE__,__func__);
				perror("Couldn't accept");
			} else {
				for(i=0;i<nsock;i++) if(socklist[i]<0) break;
				if(i==nsock) {
					if(nsock==MAX_CONN) {
						(void)shutdown(sock,SHUT_RDWR);
						(void)close(sock);
						i=-1;
						(void)fprintf(stderr,"Socket connection shutdown (too many connections)\n");
					} else {
						nsock++;
					}
				}
				if(i>=0) socklist[i]=sock;
			}
		}
		for(i=0;i<nsock;i++) if(FD_ISSET(socklist[i],&readfds)) {
			do {
				if(getppid()==1) break;
				er=0;
				j=read(socklist[i],buf1,8);
				k=0;
				while(cmdlist[k]) {
					if(!(strncmp(cmdlist[k],buf1,8))) break;
					k++;
				}
				if(j==8) switch(k) {
				 case 0: /* LMONSTAT */
					buflen=snprintf(buf,256,"%d %d %d %g %g %g %g\n",lpar->it,
										 lpar->nq,lpar->nq1,z,z1,z2,z3);
					(void)write(socklist[i],buf,buflen);
					break;
				 case 1: /* LMONINFO */
					buflen=snprintf(buf,256,"%d %d %d %d %d\n",lpar->num_iter,
										 lpar->sample_from[0],lpar->sample_from[1],
										 lpar->sample_freq[0],lpar->sample_freq[1]);
					(void)write(socklist[i],buf,buflen);
					break;
				 case 2: /* LSETINFO */
					FD_ZERO(&tmpfds);
					FD_SET(socklist[i],&tmpfds);
					nfd=select(FD_SETSIZE,&tmpfds,0,0,&timeout);
					if(nfd>0) {
						j=read(socklist[i],buf,255);
						if(j<=0) er=1;
						else {
							buf[j]=0;
							j=get_num(buf,v,5,mask);
							if(j==5) {
								if(mask[0]) lpar->num_iter=v[0];
								if(mask[1]) lpar->sample_from[0]=v[1];
								if(mask[2]) lpar->sample_from[1]=v[2];
								if(mask[3]) lpar->sample_freq[0]=v[3];
								if(mask[4]) lpar->sample_freq[1]=v[4];
							}
						}
					}
					break;
				 case 3: /* LMON_PID */
					buflen=snprintf(buf,256,"%d\n",(int)par_id);
					(void)write(socklist[i],buf,buflen);
					break;
				 case 4: /* DBR_INFO */
					if(!lpar->dbr_flag) {
						lpar->command=LMON_START_DBR;
						(void)write(socklist[i],"__END__\n",8);
					} else {
						if(n_child>=MAX_CHILD) {
							fprintf(stderr,"Too many children, can't fork()\n");
							er=1;
						} else {
							child=fork();
							if(!child) {
								for(j=0;j<nsock;j++) if(j!=i) {
									(void)shutdown(socklist[j],SHUT_RDWR);
									(void)close(socklist[j]);
								}
								(void)shutdown(s,SHUT_RDWR);
								(void)close(s);
								if(!(init_dbr_count())) write_dbr_count(socklist[i]);
								(void)write(socklist[i],"__END__\n",8);
								_exit(EXIT_SUCCESS);
							}
							er=0;
							child_id[n_child++]=child;
						}
					}
					break;
				 case 5: /* LSIGSTOP */
					(void)kill(par_id,SIGSTOP);
					break;
				 case 6: /* LSIGCONT */
					(void)kill(par_id,SIGCONT);
					break;
				 case 7: /* LSIGQUIT */
					(void)kill(par_id,SIGQUIT);
					break;
				 case 8: /* LSIGKILL */
					(void)kill(par_id,SIGKILL);
					break;
				 default:
					er=1;
 				} else er=1;
				if(er) {
					(void)shutdown(socklist[i],SHUT_RDWR);
					(void)close(socklist[i]);
					socklist[i--]=socklist[--nsock];
					break;
				}
				FD_ZERO(&tmpfds);
				FD_SET(socklist[i],&tmpfds);
				k=select(FD_SETSIZE,&tmpfds,0,0,&timeout);
			} while(k>0);
		}
	}
}

static void UnlinkFiles(void)
{
	(void)unlink(monitor_sock_addr);
	(void)unlink(loki_pid_file);
}

void start_monitor(void)
{
	int s,i,j,pid=-1,flg=0;
	struct sockaddr_un sock;
	mode_t omask;
	struct stat sbuf;
	socklen_t len;
	fd_set fds;
	struct timeval timeout;
	char buf[256],*p=0;
	char host[MAXHOSTNAMELEN+1];
	FILE *fptr;
	
#if HAVE_SYS_SYSTEMINFO_H
	if(sysinfo(SI_HOSTNAME,host,MAXHOSTNAMELEN)<0)
#else
	  if(gethostname(host,MAXHOSTNAMELEN)<0)
#endif
		 (void)strncpy(host,"UNKNOWN",MAXHOSTNAMELEN+1);
	if((s=socket(AF_UNIX,SOCK_STREAM,0))<0) {
		perror("Couldn't create socket");
		ABT_FUNC(AbMsg);
	}
	sock.sun_family=AF_UNIX;
	(void)strncpy(sock.sun_path,monitor_sock_addr,sizeof(sock.sun_path));
	omask=umask(077);
	len=sizeof(sock.sun_family)+strlen(sock.sun_path)+1;
	if(!stat(monitor_sock_addr,&sbuf)) {
		(void)fprintf(stderr,"The file '%s' is present in the current directory.\n",monitor_sock_addr);
		if(connect(s,(struct sockaddr *)&sock,len)<0) {
			if(errno==ENOTSOCK) {
				(void)fprintf(stderr,"'%s' is not a socket.\n",monitor_sock_addr);
				flg=1;
			} else {
				perror("Couldn't connect to existing loki process");
			}
		} else {
			(void)fputs("Connected to socket - trying to get process id of existing loki process\n",stderr);
			timeout.tv_sec=1;
			timeout.tv_usec=0;
			FD_ZERO(&fds);
			FD_SET(s,&fds);
			i=select(FD_SETSIZE,0,&fds,0,&timeout);
			if(i<0) perror("Couldn't write to socket");
			else if(i>0 && write(s,"LMON_PID",8)==8) {
				timeout.tv_sec=1;
				i=select(FD_SETSIZE,&fds,0,0,&timeout);
				if(i>0 && (j=read(s,buf,255))>0) {
					buf[j]=0;
					pid=atoi(buf);
				}
			}
		}
		if(pid<0 && !flg) {
			if(!(fptr=fopen(loki_pid_file,"r"))) {
				perror("Couldn't open pid file for reading");
			} else {
				if(fgets(buf,256,fptr)) {
					i=strlen(buf);
					if(buf[i-1]=='\n') buf[i-1]=0;
					i=(int)strtol(buf,&p,10);
					if(*(p++)==':') pid=i;
				}
			}
		}
		if(pid>=0) {
			(void)fprintf(stderr,"Existing loki process (%d) is%srunning on %s in this directory\n",pid,p?" apparently ":" ",p?p:host);
			(void)fprintf(stderr,"If this is not true, delete '%s' and retry\n",monitor_sock_addr);
		} else {
			if(flg==1) (void)fprintf(stderr,"Delete spurious file '%s' and retry\n",monitor_sock_addr);
			else {
				(void)fputs("Check to make sure no loki process is running locally or remotely in this directory.\n",stderr);
				(void)fprintf(stderr,"If not, delete spurious file '%s' and retry\n",monitor_sock_addr);
			}
		}
		ABT_FUNC(AbMsg);
	} else if(errno!=ENOENT) {
		perror("Problem creating socket file");
		ABT_FUNC(AbMsg);
	}
	if(atexit(UnlinkFiles)) ABT_FUNC("Unable to register exit function UnlinkFiles()\n");
	(void)umask(omask);
	if(bind(s,(struct sockaddr *)&sock,len)<0) {
		perror("Couldn't bind socket");
		ABT_FUNC(AbMsg);
	}
	if(listen(s,MAX_CONN)<0) {
		perror("Couldn't listen on socket");
		ABT_FUNC(AbMsg);
	}
	if(!(fptr=fopen(loki_pid_file,"w"))) {
		perror("Couldn't open pid file for writing");
		ABT_FUNC(AbMsg);
	}
	(void)fprintf(fptr,"%d:%s\n",(int)getpid(),host);
	(void)fclose(fptr);
	monitor_pid=fork();
	if(!monitor_pid) { /* Child */
		if(setsid()<0) {
			perror("Couldn't set new process session");
			ABT_FUNC(AbMsg);
		}
		monitor(s);
		/* Should never return */
		_exit(EXIT_FAILURE);
	} else if(monitor_pid>0) { /* Parent */
		(void)close(s);
	} else {
		perror("Couldn't fork");
		ABT_FUNC(AbMsg);
	}
	child_alive=1;
}
