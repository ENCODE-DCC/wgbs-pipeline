#ifndef _LOKI_MONITOR_H_
#define _LOKI_MONITOR_H_

void start_monitor(void);
void reaper(int);
void ignore_handler(int);

#define LMON_WIN_SIZE 64
#define LMON_MAGIC 17062000
#define LMON_START_DBR 1

#endif
