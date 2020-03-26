#ifndef POSCACHE_GUARD
#define POSCACHE_GUARD

#include "circbuf.hpp"

namespace walker {

template<class T>
struct pileup {
   uint64_t pos;
   T contents;
};

template<class T>
class pos_cache : public static_circbuf<struct pileup<T>> {
   public:
   bool contains(uint64_t pos) {
      return static_circbuf<struct pileup<T>>::at(pos).pos == pos;
   }

   void insert(uint64_t pos, T val) {
      static_circbuf<struct pileup<T>>::insert(pos, (struct pileup<T>) { pos, val });
   }

   T& at(uint64_t pos) {
      return static_circbuf<struct pileup<T>>::at(pos).contents;
   }
};

}

#endif
