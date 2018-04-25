void *_lk_malloc(size_t,const char *,char *,int);
void *_lk_calloc(size_t,size_t,const char *,char *,int);
void *_lk_realloc(void *,size_t,const char *,char *,int);

#define lk_malloc(x) _lk_malloc((x),__func__,__FILE__,__LINE__)
#define lk_calloc(x,y) _lk_calloc((x),(y),__func__,__FILE__,__LINE__)
#define lk_realloc(p,x) _lk_realloc((p),(x),__func__,__FILE__,__LINE__)
