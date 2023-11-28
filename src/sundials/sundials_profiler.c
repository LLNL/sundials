/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#include <sundials/sundials_config.h>

#if SUNDIALS_MPI_ENABLED
#include <mpi.h>
#include <sundials/sundials_mpi_types.h>
#endif

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
/* Minimum POSIX version needed for struct timespec and clock_monotonic */
#if !defined(_POSIX_C_SOURCE) || (_POSIX_C_SOURCE < 199309L)
#define _POSIX_C_SOURCE 199309L
#endif
#include <stddef.h>
#include <time.h>
#include <unistd.h>
#elif defined(WIN32) || defined(_WIN32)
#include <windows.h>
#else
#error SUNProfiler needs POSIX or Windows timers
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_types.h>

#include "sundials_debug.h"
#include "sundials_hashmap.h"

#define SUNDIALS_ROOT_TIMER ((const char*)"From profiler epoch")

#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
typedef struct timespec sunTimespec;
#else
typedef struct _sunTimespec
{
  long int tv_sec;
  long int tv_nsec;
} sunTimespec;
#endif

/* Private functions */
#if SUNDIALS_MPI_ENABLED
static int sunCollectTimers(SUNProfiler p);
#endif
static void sunPrintTimers(int idx, SUNHashMapKeyValue kv, FILE* fp, void* pvoid);
static int sunCompareTimes(const void* l, const void* r);
static int sunclock_gettime_monotonic(sunTimespec* tp);

/*
  sunTimerStruct.
  A private structure holding timing information.
 */

struct _sunTimerStruct
{
  sunTimespec* tic;
  sunTimespec* toc;
  double average;
  double maximum;
  double elapsed;
  long count;
};

typedef struct _sunTimerStruct sunTimerStruct;

static sunTimerStruct* sunTimerStructNew()
{
  sunTimerStruct* ts = (sunTimerStruct*)malloc(sizeof(sunTimerStruct));
  ts->tic            = (sunTimespec*)malloc(sizeof(sunTimespec));
  ts->toc            = (sunTimespec*)malloc(sizeof(sunTimespec));
  ts->tic->tv_sec    = 0;
  ts->tic->tv_nsec   = 0;
  ts->elapsed        = 0.0;
  ts->average        = 0.0;
  ts->maximum        = 0.0;
  ts->count          = 0;
  return ts;
}

static void sunTimerStructFree(void* TS)
{
  sunTimerStruct* ts = (sunTimerStruct*)TS;
  if (ts)
  {
    if (ts->tic) free(ts->tic);
    if (ts->toc) free(ts->toc);
    free(ts);
  }
}

static void sunStartTiming(sunTimerStruct* entry)
{
  sunclock_gettime_monotonic(entry->tic);
}

static void sunStopTiming(sunTimerStruct* entry)
{
  long s_difference  = 0;
  long ns_difference = 0;

  sunclock_gettime_monotonic(entry->toc);

  s_difference  = entry->toc->tv_sec - entry->tic->tv_sec;
  ns_difference = entry->toc->tv_nsec - entry->tic->tv_nsec;
  if (ns_difference < 0)
  {
    s_difference--;
    ns_difference = 1000000000 + entry->toc->tv_nsec - entry->tic->tv_nsec;
  }

  entry->elapsed += ((double)s_difference) + ((double)ns_difference) * 1e-9;
  entry->average = entry->elapsed;
  entry->maximum = entry->elapsed;
}

static void sunResetTiming(sunTimerStruct* entry)
{
  entry->tic->tv_sec  = 0;
  entry->tic->tv_nsec = 0;
  entry->toc->tv_sec  = 0;
  entry->toc->tv_nsec = 0;
  entry->elapsed      = 0.0;
  entry->average      = 0.0;
  entry->maximum      = 0.0;
  entry->count        = 0;
}

/*
  SUNProfiler.

  This structure holds all of the timers in a map.s
 */

struct _SUNProfiler
{
  SUNComm comm;
  char* title;
  SUNHashMap map;
  sunTimerStruct* overhead;
  double sundials_time;
};

int SUNProfiler_Create(SUNComm comm, const char* title, SUNProfiler* p)
{
  SUNProfiler profiler;
  int max_entries;
  char* max_entries_env;

  *p = profiler = (SUNProfiler)malloc(sizeof(struct _SUNProfiler));

  if (profiler == NULL) { return 0; }

  profiler->overhead = sunTimerStructNew();
  if (profiler->overhead == NULL)
  {
    free(profiler);
    *p = profiler = NULL;
    return (-1);
  }

  sunStartTiming(profiler->overhead);

  /* Check to see if max entries env variable was set, and use if it was. */
  max_entries     = 2560;
  max_entries_env = getenv("SUNPROFILER_MAX_ENTRIES");
  if (max_entries_env) max_entries = atoi(max_entries_env);
  if (max_entries <= 0) max_entries = 2560;

  /* Create the hashmap used to store the timers */
  if (SUNHashMap_New(max_entries, &profiler->map))
  {
    sunTimerStructFree((void*)profiler->overhead);
    free(profiler);
    *p = profiler = NULL;
    return (-1);
  }

  /* Attach the comm, duplicating it if MPI is used. */
#if SUNDIALS_MPI_ENABLED
  profiler->comm = SUN_COMM_NULL;
  if (comm != SUN_COMM_NULL)
  {
    MPI_Comm_dup(comm, &profiler->comm);
  }
#else
  if (comm != SUN_COMM_NULL)
  {
    free(profiler);
    return -1;
  }
  profiler->comm = SUN_COMM_NULL;
#endif

  /* Copy the title of the profiler (note strlen does not include terminating
     null character hence the +1) */
  profiler->title = malloc((strlen(title) + 1) * sizeof(char));
  strcpy(profiler->title, title);

  /* Initialize the overall timer to 0. */
  profiler->sundials_time = 0.0;

  SUNDIALS_MARK_BEGIN(profiler, SUNDIALS_ROOT_TIMER);
  sunStopTiming(profiler->overhead);

  return (0);
}

int SUNProfiler_Free(SUNProfiler* p)
{
  if (p == NULL) return (-1);

  SUNDIALS_MARK_END(*p, SUNDIALS_ROOT_TIMER);

  if (*p)
  {
    SUNHashMap_Destroy(&(*p)->map, sunTimerStructFree);
    sunTimerStructFree((void*)(*p)->overhead);
#if SUNDIALS_MPI_ENABLED
    if ((*p)->comm != SUN_COMM_NULL)
    {
      MPI_Comm_free(&(*p)->comm);
    }
#endif
    free((*p)->title);
    free(*p);
  }
  *p = NULL;

  return (0);
}

int SUNProfiler_Begin(SUNProfiler p, const char* name)
{
  int ier;
  sunTimerStruct* timer = NULL;
#ifdef SUNDIALS_DEBUG
  size_t slen;
  char* errmsg;
#endif

  if (p == NULL) return (-1);
  sunStartTiming(p->overhead);

  if (SUNHashMap_GetValue(p->map, name, (void**)&timer))
  {
    timer = sunTimerStructNew();
    ier   = SUNHashMap_Insert(p->map, name, (void*)timer);
    if (ier)
    {
#ifdef SUNDIALS_DEBUG
      slen   = strlen(name);
      errmsg = malloc(slen * sizeof(char));
      snprintf(errmsg,
               128 + slen, "(((( [ERROR] in SUNProfilerBegin: SUNHashMapInsert failed with code %d while inserting %s))))\n",
               ier, name);
      SUNDIALS_DEBUG_PRINT(errmsg);
      free(errmsg);
#endif
      sunTimerStructFree(timer);
      sunStopTiming(p->overhead);
      return (-1);
    }
  }

  timer->count++;
  sunStartTiming(timer);

  sunStopTiming(p->overhead);
  return (0);
}

int SUNProfiler_End(SUNProfiler p, const char* name)
{
  sunTimerStruct* timer;

  if (p == NULL) return (-1);
  sunStartTiming(p->overhead);

  if (SUNHashMap_GetValue(p->map, name, (void**)&timer))
  {
    sunStopTiming(p->overhead);
    return (-1);
  }

  sunStopTiming(timer);

  sunStopTiming(p->overhead);
  return (0);
}

int SUNProfiler_GetTimerResolution(SUNProfiler p, double* resolution)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  sunTimespec spec;
  clock_getres(CLOCK_MONOTONIC, &spec);
  *resolution = 1e-9 * ((double)spec.tv_nsec);

  return (0);
#elif (defined(WIN32) || defined(_WIN32))
  static LARGE_INTEGER ticks_per_sec;

  if (!ticks_per_sec.QuadPart)
  {
    QueryPerformanceFrequency(&ticks_per_sec);
    if (!ticks_per_sec.QuadPart) { return -1; }
  }

  *resolution = (double) ticks_per_sec.QuadPart;

  return (0);
#else
#error SUNProfiler needs POSIX or Windows timers
#endif
}

int SUNProfiler_GetElapsedTime(SUNProfiler p, const char* name, double* time)
{
  sunTimerStruct* timer;

  if (SUNHashMap_GetValue(p->map, name, (void**)&timer)) { return (-1); }

  *time = timer->elapsed;

  return (0);
}

int SUNProfiler_Reset(SUNProfiler p)
{
  int i                 = 0;
  sunTimerStruct* timer = NULL;

  /* Check for valid input */
  if (!p) return -1;
  if (!(p->overhead)) return -1;
  if (!(p->map)) return -1;
  if (!(p->map->buckets)) return -1;

  /* Reset the overhead timer */
  sunResetTiming(p->overhead);
  sunStartTiming(p->overhead);

  /* Reset all timers */
  for (i = 0; i < p->map->max_size; i++)
  {
    if (!(p->map->buckets[i])) continue;
    timer = p->map->buckets[i]->value;
    if (timer) sunResetTiming(timer);
  }

  /* Reset the overall timer. */
  p->sundials_time = 0.0;

  SUNDIALS_MARK_BEGIN(p, SUNDIALS_ROOT_TIMER);
  sunStopTiming(p->overhead);

  return 0;
}

int SUNProfiler_Print(SUNProfiler p, FILE* fp)
{
  int i                      = 0;
  int rank                   = 0;
  sunTimerStruct* timer      = NULL;
  SUNHashMapKeyValue* sorted = NULL;

  if (p == NULL) return (-1);
  sunStartTiming(p->overhead);

  /* Get the total SUNDIALS time up to this point */
  SUNDIALS_MARK_END(p, SUNDIALS_ROOT_TIMER);
  SUNDIALS_MARK_BEGIN(p, SUNDIALS_ROOT_TIMER);

  if (SUNHashMap_GetValue(p->map, SUNDIALS_ROOT_TIMER, (void**)&timer))
    return (-1);
  p->sundials_time = timer->elapsed;

#if SUNDIALS_MPI_ENABLED
  if (p->comm != SUN_COMM_NULL)
  {
    MPI_Comm_rank(p->comm, &rank);
    /* Find the max and average time across all ranks */
    sunCollectTimers(p);
  }
#endif

  if (rank == 0)
  {
    double resolution;
    /* Sort the timers in descending order */
    if (SUNHashMap_Sort(p->map, &sorted, sunCompareTimes)) return (-1);
    SUNProfiler_GetTimerResolution(p, &resolution);
    fprintf(fp, "\n============================================================"
                "====================================================\n");
    fprintf(fp, "SUNDIALS GIT VERSION: %s\n", SUNDIALS_GIT_VERSION);
    fprintf(fp, "SUNDIALS PROFILER: %s\n", p->title);
    fprintf(fp, "TIMER RESOLUTION: %gs\n", resolution);
    fprintf(fp, "%-40s\t %% time (inclusive) \t max/rank \t average/rank \t count \n",
            "RESULTS:");
    fprintf(fp, "=============================================================="
                "==================================================\n");

#if SUNDIALS_MPI_ENABLED
    if (p->comm == SUN_COMM_NULL)
      printf(
        "WARNING: no MPI communicator provided, times shown are for rank 0\n");
#endif

    /* Print all the other timers out */
    for (i = 0; i < p->map->size; i++)
      if (sorted[i]) sunPrintTimers(i, sorted[i], fp, (void*)p);
    free(sorted);
  }

  sunStopTiming(p->overhead);

  if (rank == 0)
  {
    /* Print out the total time and the profiler overhead */
    fprintf(fp, "%-40s\t %6.2f%% \t         %.6fs \t -- \t\t -- \n",
            "Est. profiler overhead", p->overhead->elapsed / p->sundials_time,
            p->overhead->elapsed);

    /* End of output */
    fprintf(fp, "\n");
  }

  return (0);
}

#if SUNDIALS_MPI_ENABLED
static void sunTimerStructReduceMaxAndSum(void* a, void* b, int* len,
                                          MPI_Datatype* dType)
{
  sunTimerStruct* a_ts = (sunTimerStruct*)a;
  sunTimerStruct* b_ts = (sunTimerStruct*)b;
  int i;
  for (i = 0; i < *len; ++i)
  {
    b_ts[i].average += a_ts[i].elapsed;
    b_ts[i].maximum = SUNMAX(a_ts[i].maximum, b_ts[i].maximum);
  }
}

/* Find the max and average time across all ranks */
int sunCollectTimers(SUNProfiler p)
{
  int i, rank, nranks;

  MPI_Comm comm = p->comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nranks);

  sunTimerStruct** values = NULL;

  /* Extract the elapsed times from the hash map */
  SUNHashMap_Values(p->map, (void***)&values, sizeof(sunTimerStruct));
  sunTimerStruct* reduced =
    (sunTimerStruct*)malloc(p->map->size * sizeof(sunTimerStruct));
  for (i = 0; i < p->map->size; ++i) reduced[i] = *values[i];

  /* Register MPI datatype for sunTimerStruct */
  MPI_Datatype tmp_type, MPI_sunTimerStruct;
  const int block_lens[2]     = {5, 1};
  const MPI_Datatype types[2] = {MPI_DOUBLE, MPI_LONG};
  const MPI_Aint displ[2]     = {offsetof(sunTimerStruct, tic),
                                 offsetof(sunTimerStruct, count)};
  MPI_Aint lb, extent;

  MPI_Type_create_struct(2, block_lens, displ, types, &tmp_type);
  MPI_Type_get_extent(tmp_type, &lb, &extent);
  extent = sizeof(sunTimerStruct);
  MPI_Type_create_resized(tmp_type, lb, extent, &MPI_sunTimerStruct);
  MPI_Type_commit(&MPI_sunTimerStruct);

  /* Register max and sum MPI reduction operations for our datatype */
  MPI_Op MPI_sunTimerStruct_MAXANDSUM;
  MPI_Op_create(sunTimerStructReduceMaxAndSum, 1, &MPI_sunTimerStruct_MAXANDSUM);

  /* Compute max and average time across all ranks */
  if (rank == 0)
  {
    MPI_Reduce(MPI_IN_PLACE, reduced, p->map->size, MPI_sunTimerStruct,
               MPI_sunTimerStruct_MAXANDSUM, 0, comm);
  }
  else
  {
    MPI_Reduce(reduced, reduced, p->map->size, MPI_sunTimerStruct,
               MPI_sunTimerStruct_MAXANDSUM, 0, comm);
  }

  /* Cleanup custom MPI datatype and operations */
  MPI_Type_free(&tmp_type);
  MPI_Type_free(&MPI_sunTimerStruct);
  MPI_Op_free(&MPI_sunTimerStruct_MAXANDSUM);

  /* Update the values that are in this rank's hash map. */
  for (i = 0; i < p->map->size; ++i)
  {
    values[i]->average = reduced[i].average / (double)nranks;
    values[i]->maximum = reduced[i].maximum;
  }

  free(reduced);
  free(values);

  return (0);
}
#endif

/* Print out the: timer name, percentage of exec time (based on the max),
   max across ranks, average across ranks, and the timer counter. */
void sunPrintTimers(int idx, SUNHashMapKeyValue kv, FILE* fp, void* pvoid)
{
  SUNProfiler p      = (SUNProfiler)pvoid;
  sunTimerStruct* ts = (sunTimerStruct*)kv->value;
  double maximum     = ts->maximum;
  double average     = ts->average;
  double percent = strcmp((const char*)kv->key, (const char*)SUNDIALS_ROOT_TIMER)
                     ? maximum / p->sundials_time * 100
                     : 100;
  fprintf(fp, "%-40s\t %6.2f%% \t         %.6fs \t %.6fs \t %ld\n", kv->key,
          percent, maximum, average, ts->count);
}

/* Comparator for qsort that compares key-value pairs
   based on the maximum time in the sunTimerStruct. */
int sunCompareTimes(const void* l, const void* r)
{
  double left_max;
  double right_max;

  const SUNHashMapKeyValue left  = *((SUNHashMapKeyValue*)l);
  const SUNHashMapKeyValue right = *((SUNHashMapKeyValue*)r);

  if (left == NULL && right == NULL) return (0);
  if (left == NULL) return (1);
  if (right == NULL) return (-1);

  left_max  = ((sunTimerStruct*)left->value)->maximum;
  right_max = ((sunTimerStruct*)right->value)->maximum;

  if (left_max < right_max) return (1);
  if (left_max > right_max) return (-1);
  return (0);
}

int sunclock_gettime_monotonic(sunTimespec* ts)
{
#if defined(SUNDIALS_HAVE_POSIX_TIMERS)
  return clock_gettime(CLOCK_MONOTONIC, ts);
#elif (defined(WIN32) || defined(_WIN32))
  static LARGE_INTEGER ticks_per_sec;
  LARGE_INTEGER ticks;

  if (!ticks_per_sec.QuadPart)
  {
    QueryPerformanceFrequency(&ticks_per_sec);
    if (!ticks_per_sec.QuadPart) { return -1; }
  }

  QueryPerformanceCounter(&ticks);

  /* QueryPerformanceCounter is ticks in microseconds */

  ts->tv_sec  = (long)(ticks.QuadPart / ticks_per_sec.QuadPart);
  ts->tv_nsec = (long)(((ticks.QuadPart % ticks_per_sec.QuadPart) * 1000000) /
                       ticks_per_sec.QuadPart);

  return 0;
#else 
#error SUNProfiler needs POSIX or Windows timers
#endif
}
