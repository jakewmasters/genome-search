/* A parallelized genome search program. 
 * Written by Jake Masters. 
 * 9/14/18
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

int verbose = 0;

/* Helpful constants */
#define ONE_MEGA (long)(1024 * 1024)
#define ONE_GIGA (long)(ONE_MEGA * 1024)
#define ONE_BILLION (double)1000000000.0

const int line_buffer_length = 1024;

typedef struct {
  char *sequence;				/* Entire sequence */
  char *seq_ptr;				/* Next location to store data */
  long max_length;				/* Max length allocated */
  long cur_length;				/* Current length */
} fasta_t;

/* Global variables */
char *pattern = NULL;
int num_threads = 1;
pthread_mutex_t shared_counter_mutex;
fasta_t *fasta = NULL;
long match_count = 0;
long trial_count = 0;
long chunk_size = 0;

/* Create a FASTA object; allocates memory on the heap */
fasta_t *
fasta_create(long max_length)
{
  if (verbose) {
	printf("Allocate %ld bytes\n", max_length);
  }
  fasta_t *new = malloc(sizeof(fasta_t));
  new->sequence = new->seq_ptr = malloc(max_length);
  new->max_length = max_length;
  new->cur_length = 0;
  return new;
}

/* Destroy a FASTA object; deallocates memory. */
void
fasta_destroy(fasta_t *old)
{
  free(old->sequence);
  free(old);
}

/* Read a FASTA file into a FASTA structure. Can be called multiple times and
 * will append new data to whatever is already in existing structure. For
 * example, can read multiple chromosome files into a single FASTA
 * structure. Will crater on attempts to read more data from file than will fit
 * in allocated FASTA structure. Works with both .g>SEQUEN and flat text files (but
 * prefer the zipped version to save disk space!).
 */
void
fasta_read_file(char *file_name, fasta_t *fasta)
{
  gzFile gzfp = gzopen(file_name, "rb");
  if (!gzfp) {
	fprintf(stderr, "Can't open '%s' for reading\n", file_name);
	exit(1);
  }
  printf(" LOADING %s\n", file_name);
  
  /* Reads one line at a time from the FASTA file. */
  char line_buffer[line_buffer_length];
  int lines_kept = 0;
  int lines_skipped = 0;
  while ((gzgets(gzfp, line_buffer, line_buffer_length) != NULL)) {
	char *line_ptr = line_buffer;
	if (line_buffer[0] == '>') {
	  /* Line contains text annotation; skip it. */
	  lines_skipped++;
	} else {
	  /* Valid data; copy all ACGT data from line buffer. */
	  lines_kept++;
	  while (*line_ptr != '\n') {
		*fasta->seq_ptr++ = *line_ptr++;
		fasta->cur_length++;
	  }
	}
	if (fasta->cur_length + line_buffer_length > fasta->max_length) {
	  fprintf(stderr, "Read %ld bytes; fasta buffer too small (%ld bytes)\n",
			  fasta->cur_length, fasta->max_length);
	  exit(1);
	}
  }

  if (verbose) {
	printf("%s: %d lines skipped, %d lines kept, %ld total bytes\n",
		   file_name, lines_skipped, lines_kept, fasta->cur_length);
  }
  
  gzclose(gzfp);
}

/* Return the current time. */
double
now(void)
{
  struct timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  return current_time.tv_sec + (current_time.tv_nsec / ONE_BILLION);
}

/* Return the smaller of two pointers. */
const char *
min(const char *a, const char *b)
{
  return a < b ? a : b;
}

/* Return the larger of two pointers. */
const char *
max(const char *a, const char *b)
{
  return a > b ? a : b;
}

/* Clamp 'value' to be in the range ['low' .. 'high'] */
const char *
clamp(const char *value, const char *low, const char *high)
{
  if (low >= high) {
	fprintf(stderr, "Bogus bounds [%p, %p]\n", low, high);
	exit(1);
  }

  return min(max(value, low), high);
}

/* Print 'n' padding characters to pretty-up the output. */
void
print_padding(const int n)
{
  for (int i = 0;  i < n;  i++) {
	printf(" ");
  }
}

/* Print 'length' bytes starting at 'current' from within the current data in
 * a FASTA structure. Also prints 'padding_bytes' bytes before and after the
 * range of values for context, as well as the offset into the sequence. Takes
 * care not to blow past either end of the sequence data in the FASTA
 * structure.
 */
void
bytes_around(fasta_t *fasta, char *current, int length)
{
  const int padding_bytes = 8;
  const char *fasta_last = fasta->sequence + fasta->cur_length;

  const char *first = clamp(current - padding_bytes, fasta->sequence, fasta_last);
  const char *last = clamp(current + length + padding_bytes, fasta->sequence, fasta_last);

  print_padding(padding_bytes - (current - first));

  for (const char *p = first;  p < last;  p++) {
	if (p == current) {
	  printf("[");
	}
	printf("%c", *p);
	if (p == current + length - 1) {
	  printf("]");
	}
  }

  print_padding(padding_bytes - (last - (current + length)));
  printf("%15ld\n", current - fasta->sequence);
}

/* My parallel implementation of match()
 */
void *
parallel_match(void *ptr)
{
  int pattern_length = strlen(pattern);
  char *cur_location = ptr;
  char *last_location = ptr + chunk_size-1;
  long local_count = 0;
  long local_trial = 0;

  while (cur_location <= last_location) {
	local_trial++;
	if (strncmp(cur_location, pattern, pattern_length) == 0) {
	  if (verbose) {
		bytes_around(fasta, cur_location, pattern_length);
	  }
	  local_count++;
	}
	cur_location++;
  }

  // mutex stuff! once we have the local count
  pthread_mutex_lock(&shared_counter_mutex);
  match_count += local_count;
  trial_count += local_trial;
  pthread_mutex_unlock(&shared_counter_mutex);

  return (void *)NULL;
}

/* Print a usage message and exit. */
void
usage(char *prog_name)
{
  fprintf(stderr, "%s: [-v] -b <B>|-m <MB>|-g <GB> -p <pattern> <fastafile>...\n", prog_name);
  fprintf(stderr, "  -v           enable verbose output\n");
  fprintf(stderr, "  -b <B>       allocate <B> bytes for FASTA data\n");
  fprintf(stderr, "  -m <MB>      allocate <MB> megabytes for FASTA data\n");
  fprintf(stderr, "  -g <GB>      allocate <GB> gigabytes for FASTA data\n");
  fprintf(stderr, "  -p <pattern> pattern for search [required]\n");
  fprintf(stderr, "  -h, -?       print this help and exit\n");
  fprintf(stderr, "One of -b, -m, or -g must be provided");
  fprintf(stderr, "One or more <fastafile> must appear; can be text or .gz file\n");
  fprintf(stderr, "If multiple <fastafile>s, will be concatenated and searched\n");
  exit(1);
}

void
check_thread_rtn(char *msge, int rtn) {
  if (rtn) {
    fprintf(stderr, "ERROR: %s (%d)\n", msge,rtn);
    exit(1);
  }
}

int
main(int argc, char **argv)
{
  char *prog_name = argv[0];
  long fasta_max_length = 0;

  /* Process command-line arguments; see 'man 3 getopt'. */
  int ch;
  while ((ch = getopt(argc, argv, "b:hm:g:p:vn:")) != -1) {
	switch (ch) {
	case 'b':
	  fasta_max_length = atol(optarg);
	  break;
	case 'g':
	  fasta_max_length = atol(optarg) * ONE_GIGA;
	  break;
	case 'm':
	  fasta_max_length = atol(optarg) * ONE_MEGA;
	  break;
	case 'p':
	  pattern = optarg;
	  break;
	case 'v':
	  verbose = 1;
	  break;
    case 'n':
      num_threads = atoi(optarg);
      break;
	case 'h':
	case '?':
	default:
	  usage(prog_name);
	}
  }
  argc -= optind;
  argv += optind;

  if (fasta_max_length == 0 || pattern == NULL) {
	usage(prog_name);
  }

  /* Create FASTA structure with the given length. */
  fasta = fasta_create(fasta_max_length);

  /* For each <fastafile> argument, read its data into the FASTA structure. */
  for (int idx = 0;  idx < argc;  idx++) {
	fasta_read_file(argv[idx], fasta);
  }

  // parallel stuff
  pthread_t threads[num_threads];

  int rtn = pthread_mutex_init(&shared_counter_mutex, NULL);
  check_thread_rtn("mutex init", rtn);
  
  chunk_size = fasta->cur_length / num_threads;
  void *match_ptr = fasta->sequence;

  printf("MATCHING ...\n");
  double start_time = now();
  for (int i = 0;  i < num_threads;  ++i) {
	rtn = pthread_create(&threads[i], NULL, parallel_match, match_ptr);
	check_thread_rtn("create", rtn);
    match_ptr += chunk_size;
  }

  for (int i=0; i < num_threads; ++i){
      rtn = pthread_join(threads[i], NULL);
      check_thread_rtn("join", rtn);
  }
  
  printf("    TOOK %5.3f seconds\n", now() - start_time);
  printf("   TRIED %e matches\n", (double)trial_count);
  printf(" PATTERN %s\n", pattern);
  printf("   MATCH %ld time%s\n", match_count, match_count == 1 ? "" : "s");
  
  /* Clean up and be done. */
  fasta_destroy(fasta);
  exit(0);
}
