#include "circbuf.hpp"

namespace walker {

template<class T>
typedef struct p {
   uint32_t pos;
   T contents;
} pileup_t;

template<class T>
class pos_cache : public walker::static_circbuf<pileup_t<T>> {
   public:
   bool contains(uint64_t pos) {
      return this->at(pos).pos == pos;
   }

   T& at(uint64_t pos) {
      return this->at(pos).contents;
   }
}

}
