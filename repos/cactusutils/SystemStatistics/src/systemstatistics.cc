#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <assert.h>
#include <errno.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

// There doesn't seem to be a Cactus macro defined for the Mach header
// files we actually need, so use the mach_time.h macro instead
#ifdef HAVE_MACH_MACH_TIME_H
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

#ifndef HAVE_MALLINFO

// Provide a dummy mallinfo function if none is available
struct mallinfo_t {
  int arena;
  int ordblks;
  int smblks;
  int hblks;
  int hblkhd;
  int usmblks;
  int fsmblks;
  int uordblks;
  int fordblks;
  int keepcost;
};

struct mallinfo_t mallinfo()
{
  struct mallinfo_t m;
  m.arena = 0;
  m.ordblks = 0;
  m.smblks = 0;
  m.hblks = 0;
  m.hblkhd = 0;
  m.usmblks = 0;
  m.fsmblks = 0;
  m.uordblks = 0;
  m.fordblks = 0;
  m.keepcost = 0;
  return m;
}

#endif

#ifndef _MACH_INIT_

static unsigned long int get_rss()
{
  unsigned int size=0; //       total program size
  unsigned int resident=0;//   resident set size
  unsigned int share=0;//      shared pages
  unsigned int text=0;//       text (code)
  unsigned int lib=0;//        library
  unsigned int data=0;//       data/stack

  int page_size = sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  // If the /proc filesystem does not exist, this file will not be
  // found and the function will return 0
  if (pf) 
  {
    if (fscanf(pf, "%u %u %u %u %u %u",
               &size, &resident, &share, &text, &lib, &data) != 6)
    {
      CCTK_WARN(1, "Error while reading memory statistics (rss); results will be invalid");
      fclose(pf);
      return 0;
    }
    
    fclose(pf);
  }
  return (unsigned long int ) resident * (unsigned long int) page_size;
}

#else

// The code to get the RSS from Mac OS has been modified from
//   http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html

static unsigned long int get_rss()
{
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    task_t task = MACH_PORT_NULL;

    if (task_for_pid(current_task(), getpid(), &task) != KERN_SUCCESS)
        abort();

    task_info(task, TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    return t_info.resident_size;
}

#endif

static unsigned int get_majflt()
{
  int pid;
  char exe[256];
  char  state;
  int dummyi;
  unsigned long int dummyu;
  unsigned long int majflt;

  unsigned int page_size = sysconf(_SC_PAGESIZE);

  char buf[30];
  snprintf(buf, 30, "/proc/%u/stat", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  // If the /proc filesystem does not exist, this file will not be
  // found and the function will return 0
  if (pf) 
  {
    if (fscanf(pf, "%d %s %c %d %d %d %d %d %lu %lu %lu %lu", 
           &pid, exe, &state, &dummyi,  &dummyi, &dummyi, &dummyi, 
               &dummyi, &dummyu, &dummyu, &dummyu, &majflt) != 12)
    {
      CCTK_WARN(1, "Error while reading memory statistics (majflt); results will be invalid");
      fclose(pf);
      return 0;
    }

    fclose(pf);
  }
  return majflt * page_size;
}

static long long int get_swap_kB()
{
  FILE *f = fopen("/proc/meminfo", "r");
  if (f == 0)
    return -1;
  const int buf_len = 100;
  char buffer[buf_len];
  long long int swap_total = 0;
  long long int swap_free = 0;
  bool read_swap_total = false;
  bool read_swap_free = false;

  while(!feof(f) && fgets(buffer, buf_len, f) != NULL)
  {
    char key[100];
    char unit[100];
    long long int val = -1;
    sscanf(buffer, "%s %lld %s", key, &val, unit);
    if (strcmp(key, "SwapTotal:") == 0)
    {
      assert(strcmp(unit, "kB") == 0);
      swap_total = val;
      read_swap_total = true;
    }
    else if (strcmp(key, "SwapFree:") == 0)
    {
      assert(strcmp(unit, "kB") == 0);
      swap_free = val;
      read_swap_free = true;
    }
  }

  if (!read_swap_total || !read_swap_free)
  {
    CCTK_WARN(1, "Unable to read swap usage from /proc/meminfo");
    swap_total = 0; swap_free = 0;
  }

  fclose(f);
  return swap_total - swap_free;
}

#ifdef HAVE_MALLOC_INFO
// numbers compared to mallinfo seem to be off by the size of the xml string
// rounded to full pages plus one page, plus 16 byte. No idea where the 16 byte
// come from.
#define FUDGE_AMOUNT 16

static size_t extract_value(const char* buf, const char* fmt)
{
  size_t value = 0;
  bool found = false;
  /* (p = strchr(p, '\n')) && ++p avoids SEGFAULT if '\n' is not found */
  for(const char* p = buf ; p && *p ; (p = strchr(p, '\n')) && ++p) {
    /* skip over all individual heaps until we find the true totals */
    int nr;
    if(sscanf(p, "<heap nr=\"%d\">", &nr)) {
      for(/*nop*/ ; p && *p ; (p = strchr(p, '\n')) && ++p) {
        const char* endheap = "</heap>";
        if(strncmp(p, endheap, strlen(endheap)) == 0) {
          break;
        }
      }
      /* reached end of string without seeing a '\n', this shouldn't really happen */
      if(!p || !*p)
        break;
      /* current line is </heap> which never matches */
      continue;
    }

    found = sscanf(p, fmt, &value);
    if(found)
      break;
  }
  if(!found) {
    CCTK_VWARN(CCTK_WARN_COMPLAIN, "Could not extract using format '%s' from '%s'",
               fmt, buf);
  }
  return value;
}

static char* get_malloc_info(size_t* sz)
{
  char* buf = NULL;
  FILE* memfh = open_memstream(&buf, sz);
  if(memfh != NULL) {
    const int ierr = malloc_info(0, memfh);
    fclose(memfh); // must be before call to free()
    if(ierr != 0) {
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "malloc_info failed: %s", strerror(ierr));
      free(buf);
      buf = NULL;
    }

    const size_t malloc_version = extract_value(buf, "<malloc version=\"%zu\"/>");
    static bool have_warned = false;
    if(!have_warned && malloc_version != 0 && malloc_version != 1) {
      CCTK_VWARN(CCTK_WARN_COMPLAIN, "Unexpected malloc version: %zu, only know how to handle 1",
                 malloc_version);
      have_warned = true;
    }
  } else {
    CCTK_VWARN(CCTK_WARN_COMPLAIN, "Could not open memory file handle: %s",
               strerror(errno));
  }
  return buf;
}

static size_t get_uordblks()
{
  size_t uordblks = 0;
  const long page_size = sysconf(_SC_PAGE_SIZE);

  size_t sz; /* keep track of how much memory is used by us */
  char* buf = get_malloc_info(&sz);
  if(buf != NULL) {
    const size_t total_aspace = extract_value(buf, "<aspace type=\"total\" size=\"%zu\"/>");
    const size_t total_fastavail =
      extract_value(buf, "<total type=\"fast\" count=\"%*zu\" size=\"%zu\"/>");
    const size_t total_avail =
      extract_value(buf, "<total type=\"rest\" count=\"%*zu\" size=\"%zu\"/>");
    // Fudge free amount a bit to account for memory used by XML string buffer
    const size_t fudge = ((sz + page_size-1) & ~(page_size-1)) + page_size + FUDGE_AMOUNT;
    uordblks = total_aspace - total_avail - total_fastavail - fudge;
  }
  free(buf);

  return uordblks;
}

static size_t get_hblkhd()
{
  size_t hblkhd = 0;

  size_t sz; /* keep track of how much memory is used by us */
  char* buf = get_malloc_info(&sz);
  if(buf != NULL) {
    const size_t mmapped_mem =
      extract_value(buf, "<total type=\"mmap\" count=\"%*zu\" size=\"%zu\"/>");
    hblkhd = mmapped_mem;
  }
  free(buf);

  return hblkhd;
}
#else
static size_t get_uordblks()
{
  return mallinfo().uordblks;
}

static size_t get_hblkhd()
{
  return mallinfo().hblkhd;
}
#endif

#ifdef TEST_MALLOC_INFO
#define CCTK_VWARN(level, fmt, ...) printf(fmt, __VA_ARGS__), putchar('\n')

int main() {
  size_t total = 0;

  // try to preallocate all internal structures
  for(int j = 1 ; j < 3 ; j++) {
    size_t sz = 1;
    for(size_t i = 0 ; i < 20 ; i++) {
      malloc(sz);
      total += sz;
      sz *= 2;
    }
  }
  mallinfo();
  get_uordblks();

  #define DIM(a) (int)(sizeof(a)/sizeof(a[0]))
  // some trial sizes
  size_t sz1[] = {4096, 4096, 12, 100000};
  void* p1[DIM(sz1)];
  bool free1[DIM(sz1)] = {true,false,false,true};

  size_t sz2[] = {4096, 128, 1};
  void* p2[DIM(sz2)];
  bool free2[DIM(sz2)] = {false,false,false};

  int ref1 = mallinfo().uordblks;
  size_t mine1 = get_uordblks();
  for(int i = 0 ; i < DIM(sz1) ; i++) {
    p1[i] = malloc(sz1[i]);
    total += sz1[i];
  }
  for(int i = 0 ; i < DIM(sz1) ; i++) {
    if(free1[i]) free(p1[i]);
    total -= sz1[i];
  }
  int ref2 = mallinfo().uordblks;
  size_t mine2 = get_uordblks();
  for(int i = 0 ; i < DIM(sz2) ; i++) {
    p2[i] = malloc(sz2[i]);
    total += sz2[i];
  }
  for(int i = 0 ; i < DIM(sz2) ; i++) {
    if(free2[i]) free(p2[i]);
    total -= sz2[i];
  }
  int ref3 = mallinfo().uordblks;
  size_t mine3 = get_uordblks();

  printf("ref1: %d ref2: %d ref3: %d\n", ref1, ref2, ref3);
  printf("mine1: %zu mine2: %zu mine3: %zu\n", mine1, mine2, mine3);

  printf("mine21: %zu mine32: %zu mine31: %zu\n", mine2-mine1, mine3-mine2, mine3-mine1);
  printf("ref21: %d ref32: %d ref31: %d\n", ref2-ref1, ref3-ref2, ref3-ref1);
  printf("total: %zu\n", total);

  malloc_stats();

  return 0;
}
#endif

extern "C" void SystemStatistics_Collect(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_SystemStatistics_Collect
  DECLARE_CCTK_PARAMETERS

  const int mb = 1024*1024;
  const int kb = 1024;

  *maxrss = get_rss();
  *majflt = get_majflt();
  *arena = mallinfo().arena;
  *ordblks = mallinfo().ordblks;
  *hblks = mallinfo().hblks;
  *hblkhd = get_hblkhd();
  *uordblks = get_uordblks();
  *fordblks = mallinfo().fordblks;
  *keepcost = mallinfo().keepcost;
  *swap_used = get_swap_kB() * kb;

  *maxrss_mb = get_rss() / mb;
  *majflt_mb = *majflt / mb;
  *arena_mb = *arena / mb;
  *ordblks_mb = *ordblks / mb;
  *hblks_mb = *hblks / mb;
  *hblkhd_mb = *hblkhd / mb;
  *uordblks_mb = *uordblks / mb;
  *fordblks_mb = *fordblks / mb;
  *keepcost_mb = *keepcost / mb;
  *swap_used_mb = *swap_used / mb;

  *maxrss_kb = get_rss() / kb;
  *majflt_kb = *majflt / kb;
  *arena_kb = *arena / kb;
  *ordblks_kb = *ordblks / kb;
  *hblks_kb = *hblks / kb;
  *hblkhd_kb = *hblkhd / kb;
  *uordblks_kb = *uordblks / kb;
  *fordblks_kb = *fordblks / kb;
  *keepcost_kb = *keepcost / kb;
  *swap_used_kb = *swap_used / kb;

}
