#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
#include "TestArray.h"

static const int LABEL_SIZE = 116;

int TestArray_write_helper(
    const char *filename,
    const char *label,
    const void *array,
    int elemsize,
    int length,
    int isfp // is floating point?
    ) {
  const char *stype = "unknown";
  if      (  isfp && elemsize==8 ) stype = "double";
  else if (  isfp && elemsize==4 ) stype = "float";
  else if ( !isfp && elemsize==8 ) stype = "int64";
  else if ( !isfp && elemsize==4 ) stype = "int32";
  else if ( !isfp && elemsize==2 ) stype = "int16";
  else if ( !isfp && elemsize==1 ) stype = "int8";
  printf("Writing %s array of length %d to binary file %s\n",
      stype, length, filename);
  FILE *fp = fopen(filename, "wb");
  if ( fp == 0 ) return -1;
  char padbuf[LABEL_SIZE] = { 0 };
  strncpy(padbuf, label, LABEL_SIZE-1);
  fwrite(padbuf, 1, LABEL_SIZE, fp);
  fwrite(&isfp, sizeof(int), 1, fp);
  fwrite(&elemsize, sizeof(int), 1, fp);
  fwrite(&length, sizeof(int), 1, fp);
  fwrite(array, elemsize, length, fp);
  fclose(fp);
  fflush(stdout);
  return 0;
}

#ifdef TEST_ARRAY_STANDALONE

typedef char int8;
typedef short int16;
typedef int int32;
typedef long long int int64;

//
// Pass unitialized pointers for label and array buffer.
// Storage allocated based on file content.
//
int TestArray_read(
    const char *filename,
    char **p_label,
    void **p_buffer,
    int *p_isfp,  // set to nonzero if file indicates floating point
    int *p_elemsize,
    int *p_length
    ) {
  FILE *fp = fopen(filename, "rb");
  if ( fp == 0 ) return -1;
  char *label = (char *) malloc(LABEL_SIZE);
  int isfp, elemsize, length;
  fread(label, 1, LABEL_SIZE, fp);
  fread(&isfp, sizeof(int), 1, fp);
  fread(&elemsize, sizeof(int), 1, fp);
  fread(&length, sizeof(int), 1, fp);
  void *buffer = malloc(elemsize * length);
  fread(buffer, elemsize, length, fp);
  fclose(fp);
  *p_label = label;
  *p_buffer = buffer;
  *p_isfp = isfp;
  *p_elemsize = elemsize;
  *p_length = length;
  return 0;
}

//
// Free label and array buffer memory allocated by file reader.
//
void TestArray_free_memory(
    char *label,
    void *buffer
    ) {
  free(label);
  free(buffer);
}

//
// Read ith of sz-sized elements from buf.  Return as double.
//
double TestArray_read_elem_double(void *buf, int i, int sz) {
  double x;
  void *p = (void *) ((const char *) buf + i*sz);
  switch (sz) {
    case 4:
      x = *((float *) p);
      break;
    case 8:
      x = *((double *) p);
      break;
    default:
      fprintf(stderr,
          "Unrecognized element size %d for floating point array\n", sz);
      exit(1);
  }
  return x;
}

//
// Read ith of sz-sized elements from buf.  Return as double.
//
int64 TestArray_read_elem_int(void *buf, int i, int sz) {
  int64 n;
  void *p = (void *) ((const char *) buf + i*sz);
  switch (sz) {
    case 1:
      n = *((int8 *) p);
      break;
    case 2:
      n = *((int16 *) p);
      break;
    case 4:
      n = *((int32 *) p);
      break;
    case 8:
      n = *((int64 *) p);
      break;
    default:
      fprintf(stderr,
          "Unrecognized element size %d for integer array\n", sz);
      exit(1);
  }
  return n;
}

int main(int argc, const char *argv[]) {
  if (sizeof(int8) != 1 ||
      sizeof(int16) != 2 ||
      sizeof(int32) != 4 ||
      sizeof(int64) != 8) {
    fprintf(stderr, "Integer sizes are incorrect. Must rebuild.\n");
    exit(1);
  }
  if (argc < 3) {
    fprintf(stderr, "Compare two binary files created using TestArray.\n");
    fprintf(stderr, "Syntax: %s file1 file2 [first_index] [last_index]\n",
      argv[0]);
    exit(1);
  }
  const char *fname1 = argv[1];
  const char *fname2 = argv[2];
  int first = 0;
  int last = INT_MAX;
  if (argc > 3) first = atoi(argv[3]);
  if (argc > 4) last = atoi(argv[4]);

  char *label1, *label2;
  void *buffer1, *buffer2;
  int isfp1, isfp2;
  int elemsize1, elemsize2;
  int length1, length2;

  if ( TestArray_read(fname1, &label1, &buffer1,
        &isfp1, &elemsize1, &length1) ) {
    fprintf(stderr, "Unable to read first file %s\n", fname1);
    exit(1);
  }

  if ( TestArray_read(fname2, &label2, &buffer2,
        &isfp2, &elemsize2, &length2) ) {
    fprintf(stderr, "Unable to read second file %s\n", fname2);
    exit(1);
  }

  if ( isfp1 != isfp2 ) {
    fprintf(stderr, "Trying to compare floating point and integer data\n");
    exit(1);
  }

  if (last >= length1) last = length1-1;
  if (last >= length2) last = length2-1;

  printf("%s\n%s\n", label1, label2);
  if (isfp1) {
    printf("%6s  %14s  %14s  %14s  %14s\n",
        "index", "file1", "file2", "abserr", "relerr");
    for (int i=first;  i <= last;  i++) {
      double a = TestArray_read_elem_double(buffer1, i, elemsize1);
      double b = TestArray_read_elem_double(buffer2, i, elemsize2);
      double abserr = a - b;
      double relerr = (a != 0 ? fabs(abserr) / fabs(a) : 0);
      printf("%6d  %14.8f  %14.8f  %14.8f  %14e\n",
          i, a, b, abserr, relerr);
    }
  }
  else {
    printf("%6s  %14s  %14s  %14s\n",
        "index", "file1", "file2", "abserr");
    for (int i=first;  i <= last;  i++) {
      int64 a = TestArray_read_elem_int(buffer1, i, elemsize1);
      int64 b = TestArray_read_elem_int(buffer2, i, elemsize2);
      int64 abserr = a - b;
      printf("%6d  %14ld  %14ld  %14ld\n", i, a, b, abserr);
    }
  }

  TestArray_free_memory(label1, buffer1);
  TestArray_free_memory(label2, buffer2);

  return 0;
}

#endif // TEST_ARRAY_STANDALONE
