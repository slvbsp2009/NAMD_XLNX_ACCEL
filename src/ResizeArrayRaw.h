/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
  Data storage class for ResizeArray and related container classes.
  Copy and assignment copy raw pointers, use copy() to copy data.
  Destructor does not delete data, call free first.
*/


#ifndef RESIZEARRAYRAW_H
#define RESIZEARRAYRAW_H

#include <new>
#include <string.h>
#include "common.h"

#define ResizeArrayGrowthFactor 1.5
#define ResizeArrayMinSize 8

// Need this juju to use templated friend below
template <class Type> class ResizeArray;
template <class Type> class SortableResizeArray;
template <class Type> class SortedArray;
template <class Type> class UniqueSortedArray;
template <class Type> class ResizeArrayIter;

// Class assumes that one can bit move objects
// around on array.  This will be true
// as long as object has no pointers to itself.
template <class Elem> class ResizeArrayRaw {

  private:
    Elem *array;
    unsigned char *varray;

    // Use of int (int32) is safe for counting elements.
    int arraySize;
    int allocSize;

    // No constructor run on new elements
    // arraySize is not adjusted, only allocSize
    void resizeRaw(int size) {
      if (size <= allocSize) return;
  
      // Note that ResizeArrayGrowthFactor is type double so the
      // multiplication below is also done in double.
      // Assuming int32, the conversion from double overflows at
      //   allocSize = 1431655765
      // Conditional remains safe if overflow causes RHS < 0.
      if (size < (int)(allocSize*ResizeArrayGrowthFactor)) {
        size = (int)(allocSize*ResizeArrayGrowthFactor);
      }
      // The next conditional will test true either for very small or
      // unallocated arrays or if the above conditional causes overflow.
      if ( (size-allocSize) < ResizeArrayMinSize) {
        // Overflow occurs here only when allocSize is within
        // ResizeArrayMinSize of 2^31.
        size = allocSize+ResizeArrayMinSize;
      }

      // Align everything to 64-byte boundaries (if possible).
      // Use of sizeof below promotes expression to type size_t.
      unsigned char *tmpv = new unsigned char[size*sizeof(Elem)+63];
      //Elem *tmpa = (Elem *)((((long)tmpv)+63L)&(-64L));
      // Someday we might need this alternate form.
      //
      // The following pointer manipulation is dangerous. We should use
      // reinterpret_cast<uintptr_t> or equivalent to convert tmpv+63
      // to an unsigned integer type before doing bitwise and operation.
      //
      // The 63 value (64-byte alignment) aligns for AVX-512 registers.
      // Would be better to use a macro or maybe even a template parameter
      // instead of hard-coding the alignment.  With a template parameter,
      // we could override the alignment where optimization warrants or
      // target builds towards a particular alignment, e.g., AVX needs
      // only 32-byte alignment.
      Elem *tmpa = (Elem *)(tmpv+63 - (((long)(tmpv+63))&(63L)));
      if (arraySize) CmiMemcpy((void *)tmpa, (void *)array, sizeof(Elem)*arraySize);
  
      if (allocSize) delete[] varray;
      varray = tmpv;
      array = tmpa;
      allocSize = size;
    }

    friend class ResizeArray<Elem>;
    friend class SortableResizeArray<Elem>;
    friend class SortedArray<Elem>;
    friend class UniqueSortedArray<Elem>;
    friend class ResizeArrayIter<Elem>;

    inline int size(void) const { return arraySize; }
    inline Elem &operator[](int index) const { return array[index]; }

    // DMK - MIC Support - Allow us to see the buffer's size, not just how much of it is used
    #if NAMD_MIC != 0
      inline int bufSize(void) const { return allocSize; }
    #endif

    // Default constructor 
    ResizeArrayRaw(void) : 
      array((Elem *)0), varray((unsigned char *)0), arraySize(0), allocSize(0) { }

    // Encap a pre-existing array
    ResizeArrayRaw( Elem * * const array, int arraySize, int allocSize) {
      if (allocSize < arraySize) allocSize = arraySize;
      this->allocSize = allocSize;
      this->arraySize = arraySize;
      varray = (unsigned char *)*array;
      this->array = (Elem *)*array;
      *array = 0;
    }
  
    // copy data
    void copy(const ResizeArrayRaw<Elem> &rar ) {
  
      // Clean up this array
      resize(0);
      resizeRaw(rar.size());
  
      CmiMemcpy((void*)array, (void*)rar.array, sizeof(Elem)*rar.size());
      arraySize = rar.size();
    }
  
    // Properly constructs default object on new elements
    // Properly destructs on removed elements
    // arraySize is properly updated
    void resize(int size) {
      int i;
  
      if (size < arraySize) {
        for (i=size; i<arraySize; i++) {
          array[i].~Elem();
        }
      } else if (size > arraySize) {
        resizeRaw(size);
        for (i=arraySize; i<size; i++) {
          new ((void *)&array[i]) Elem;
        }
      }
      arraySize = size;
    }
  
    // needed for ResizeArray destructor
    void free(void) {
      for (int i=0; i < arraySize; i++) {
        array[i].~Elem();
      }
      delete[] varray;
    }
  
    // resize to 0 and free storage
    void clear(void) {
      free();
      array = 0;
      varray = 0;
      arraySize = 0;
      allocSize = 0;
    }

    inline void del(int index, int number) {
      if (number) {
        // Destruct objects to be deleted
        for (int i=index; i < index+number; i++) {
          array[i].~Elem();
        }
  
        // Shift down
        memmove((void *)(array+index),
           (void *)(array+index+number),
           (arraySize-number-index)*sizeof(Elem));
      
        // fixup size of array
        arraySize -= number;
      }
    }
  
    
    // Insert element in array
    // If index is over the end of array, default
    // constructor should be run over blank elements to pad.
    inline void ins(const Elem &e, int index) {
      // Size array depending if index is in current array or reaches beyond.
      if (index < arraySize) {
        resizeRaw(arraySize+1);
        // Shift up
        memmove((void *)(array+index+1),
          (void *)(array+index),
          (arraySize-index)*sizeof(Elem));
      } else {
        resizeRaw(index+1);
      }
      
      // Write in new element via assignment - allows any refcounting
      // etc. to take place correctly!
      new((void *)&array[index]) Elem;
      array[index] = e;
    
      // Take care of fill and setting correct arraySize 
      if (index > arraySize) {
        for (Elem *tmp = array+arraySize; tmp < array+index; tmp++) {
          new ((void *)tmp) Elem;
        }
        arraySize = index+1;
      } else {
        arraySize++;
      }
    }

    inline int find(const Elem &e) const {
      for (int i=0; i<arraySize; i++) {
        if (array[i] == e) return i;
      }
      return -1;
    }
};	// end template definition

#endif
