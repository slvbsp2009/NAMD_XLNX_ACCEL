#ifndef TEST_ARRAY_H
#define TEST_ARRAY_H

#include <type_traits>

void NAMD_die(const char *);

int TestArray_write_helper(
    const char *filename,
    const char *label,
    const void *array,
    int elemsize,
    int length,
    int isfp
    );

/**
 * Write an array to file for later analysis.
 * The array type needs to be a standard arithmetic type,
 * as designated true by std::is_arithmetic<T>.
 *
 * File format is:
 * - label is padded to 116 characters including null terminator
 * - nonzero 4-byte int if array is floating point
 * - number of bytes per element stored as 4-byte int
 * - length of the array stored as 4-byte int
 * - the array data
 *
 * Returns zero for success or nonzero for error.
 */
template <typename T>
int TestArray_write(
    const char *filename,  ///< save to file name
    const char *label,     ///< label for array (truncated to 120 characters)
    const T *array,        ///< points to the array to save
    int length             ///< length of the array
    ) {
  if ( ! std::is_arithmetic<T>::value ) {
    NAMD_die("Called TestArray_write with non-aritmetic type");
  }
  int isfp = std::is_floating_point<T>::value;
  return TestArray_write_helper(filename, label,
      (const void *)array, sizeof(T), length, isfp);
}

#endif
