#ifndef _LK_LONG_H_
#define _LK_LONG_H_

#ifdef USE_LONGLONG
typedef long long lk_long;
typedef unsigned long long lk_ulong;
#define LK_LONG_BIT 64
#define LK_LONG_MAX 0x7fffffffffffffff
#else
typedef long lk_long;
typedef unsigned long lk_ulong;
#define LK_LONG_BIT LONG_BIT
#define LK_LONG_MAX LONG_MAX
#endif

#endif
