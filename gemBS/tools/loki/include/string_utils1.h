#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

#include "utils.h"

typedef struct {
  char *s;
  int len,slen;
} string;

string *add_to_string(string *,char);
string *make_string(int,int);
string *addn_to_string(string *,char *,int);
string *paddn_to_string(string *,char *,int);
string *add_strings(string *,string *);
char *copy_cstring(string *);
string *copy_string(string *);
string *strip_quotes(string *);
string *copy_string1(string *,string *);
void free_string(string *);
string *int_to_string(int);
string *real_to_string(double);
int string_cmp(string *,string *);

#define LOG_TEN 2.30258509299404568402

#endif
