#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>

#include "xml.h"

static char *infile;

static void print_att_list(const args *arg)
{
	while(arg) {
		printf(" %s='%s'",arg->name,arg->att);
		arg=arg->next;
	}
}

static void start(const string *name,const args *attr)
{
	printf("<%s",name->s);
	if(attr) print_att_list(attr);
	printf(">");
}

static void end(const string *name)
{
	printf("</%s>",name->s);
}

static void content(const string *s)
{
	fputs(s->s,stdout);
}

static void xml_error(const int line, const char *fmt, ...)
{
	va_list a;
	
	if(stdout) (void)fflush(stdout);
	(void)fprintf(stderr,"\n[%s:%d] ",infile,line);
	va_start(a,fmt);
	(void)vfprintf(stderr,fmt,a);
	va_end(a);
}

int main(int argc,char *argv[])
{
	XML_handler handle;
	FILE *fptr;
	
	handle.start_element=start;
	handle.end_element=end;
	handle.content=content;
	handle.comment=0;
	handle.error=xml_error;
	if(argc>1) {
		infile=argv[1];
		if(!(fptr=fopen(argv[1],"r"))) exit(-1);
		Read_XML(fptr,&handle);
		fclose(fptr);
	}
	return 0;
}
