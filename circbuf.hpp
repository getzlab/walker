#ifndef CIRCBUF_GUARD
#define CIRCBUF_GUARD

#include <cstdint>
#include <cstdlib>

/*
x[pos % 1000] = 64 bit int

gpos (32 bits) [could be 28 bits for humans -- log2(chr1 length) = 27.9]
event (8 bits)
16 [18] bits for alt and refcounts

Implement with failure on collision (if buffer is big enough, this should never happen)
If we wanted to get fancy, we could make it decay to a pointer to a linked list for resolution
*/

namespace walker {

// TODO: make these all operators. make virtual where useful.
//       collision checking?
template<class T = uint64_t, uint8_t POW = 10>
class static_circbuf {
   public:
   T& at(uint64_t pos) {
      return buffer[pos & mask];
   }

   void insert(uint64_t pos, T val) {
      buffer[pos & mask] = val;
   }

   void erase(uint64_t pos) {
      buffer[pos & mask] = T();
   }

   static_circbuf() : buffer() {}

   //virtual bool equals(uint64_t pos, T val); // TODO: make this an operator

   protected:
   uint64_t size {1 << POW};
   uint64_t mask {(1 << POW) - 1};
   T buffer[1 << POW];
};

}

#endif
