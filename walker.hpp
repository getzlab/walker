#ifndef WALKER_GUARD
#define WALKER_GUARD

#include <SeqLib/BamReader.h>
#include <SeqLib/BamRecord.h>
#include <SeqLib/RefGenome.h>
#include <SeqLib/GenomicRegion.h>
#include <SeqLib/GenomicRegionCollection.h>

#include <string>
#include <vector>

namespace walker {

// taken from BWA's bntseq.c (Heng Li)
static const uint8_t pack_2bit[256] = {
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  /* A C G */
   4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  /* T */
   4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  /* a c g */
   4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  /* t */
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  // XXX: is going beyond "t" really necessary?
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// convert 4 bit (nybble) representation of base in read to 2 bit (assumes we
// don't care about N's in reads)
static const uint8_t nyb_to_2bit[16] = {4, 0, 1, 4,   /* A C */
					2, 4, 4, 4,   /* G */
					3, 4, 4, 4,   /* T */
					4, 4, 4, 4};  /* N */

using namespace std;

class walker {
   public:
   /** Check if read has edit distance zero (excluding clipped bases)
    *  @return Return false if edit distance > 0; true if nonexistent NM tag
    *          (to be conservative)
    */
   static bool EDz(const SeqLib::BamRecord& record);

   /* Get all genomic positions in read not matching reference
    */
   virtual vector<uint64_t> nonref_pos(const SeqLib::BamRecord& record);

   /* iterators */
   virtual void walk() { return; } // start-to-end
   virtual void walk(SeqLib::GenomicRegion) { return; } // single region
   virtual void walk(SeqLib::GenomicRegionCollection<>) { return; } // multiple region

   /* constructors */
   walker(const string& bam_in, const string& ref_fa);

   // members

   protected:
   SeqLib::BamReader reader;
   SeqLib::BamHeader header;
   SeqLib::BamRecord cur_read;
   //multimap<uint32_t, r_read> read_ends;
   //set<uint32_t> nonref_poses;
   SeqLib::RefGenome reference;
   //string output_stem;
   //uint32_t num_pileups_dumped;
};

}

#endif
