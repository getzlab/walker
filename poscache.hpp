#include "circbuf.hpp"

namespace walker {

template<class T>
struct pileup {
   uint32_t pos;
   T contents;
};

template<class T>
class pos_cache : public static_circbuf<struct pileup<T>> {
   public:
   bool contains(uint64_t pos) {
      return this->at(pos).pos == pos;
   }

   T& at(uint64_t pos) {
      return this->at(pos).contents;
   }
};

}
