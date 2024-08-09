#ifndef WALKER_GUARD
#define WALKER_GUARD

#include <SeqLib/BamReader.h>
#include <SeqLib/BamRecord.h>
#include <SeqLib/RefGenome.h>
#include <SeqLib/GenomicRegion.h>
#include <SeqLib/GenomicRegionCollection.h>

#include <string>
#include <vector>
#include <chrono>

namespace walker {

using namespace std;

// taken from BWA's bntseq.c (Heng Li)
static const uint8_t pack_2bit[256] = {
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
   4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,  /* - */
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

typedef struct nonref_pos_t {
   uint32_t refpos;
   uint16_t readpos;
   uint16_t jointpos;
   uint16_t cig;
} nonref_pos_t;

class walker {
   public:
   /** Check if read has edit distance zero (excluding clipped bases)
    *  @return Return false if edit distance > 0; true if nonexistent NM tag
    *          (to be conservative)
    */
   static bool EDz(const SeqLib::BamRecord& record);

   /* Get all genomic positions in read not matching reference
     * @return lower 32 bits contain absolute reference position of mismatch
     *         high 32 bits contain: 
     *           - lower 16: position along the read of mismatch
     *           - upper 16: 2 bits for CIGAR op (as defined in sam.h), 14 bits
     *             for op. length (clamped at 2^14 = 16384)
    */
   virtual vector<nonref_pos_t> nonref_pos(const SeqLib::BamRecord& record);

   /* iterators */

   // 1. we're interested in all positions
   virtual void walk(); // start-to-end
   virtual void walk(const SeqLib::GenomicRegion&); // single region
   virtual void walk(const SeqLib::GenomicRegionCollection<>&); // multiple region

   // 2. we're only interested in positions intersecting a position list
   virtual void walk(uint8_t*) { } // start-to-end
   virtual void walk(const SeqLib::GenomicRegion&, uint8_t*) { } // single region
   virtual void walk(const SeqLib::GenomicRegionCollection<>&, uint8_t*) { } // multiple region

   /** Function to apply to each read in iterator
    *  @return Return false if iterator should break; true otherwise
    */
   virtual bool walk_apply(const SeqLib::BamRecord& record) { return false; };

   /** Increment current position, accounting for chromosome boundaries
    *  @parameter curchr: current chromosome index (with respect to header order)
    *  @parameter curpos: position within chromosome
    */
   void increment_pos(uint16_t& curchr, uint32_t& curpos);

   /** Check if read has edit distance zero (excluding clipped bases)
    *  @return Return true if read fails any filters; false otherwise
    */
   virtual bool filter_read(const SeqLib::BamRecord& record);

   /** Print current status:
    *  - Current position of this.reader
    *  - Number of reads per second consumed since print_status() was last called
    */
   void print_status();

   /** Set current output file
    *  @parameter String pointing to output file (or "-" for stdout)
    */
   bool set_output_file(const string& outfile);

   /** Close current output file. Does not reset this->outfile_name.
    */
   bool close_output_file();

   /* constructors */
   walker(const string& bam_in, const string& ref_fa = "");
   ~walker();

   // members

   protected:
   // SeqLib objects
   SeqLib::BamReader reader;
   SeqLib::BamHeader header;
   SeqLib::BamRecord cur_read;
   SeqLib::RefGenome reference;

   // output file
   FILE* outfile = NULL;
   string outfile_name;

   // total number of reads consumed (now, at previous query)
   uint64_t n_reads = 0, n_reads_last = 0;

   // total number of reads actually processed
   // must be incremented by user in walker::walk()
   uint64_t n_reads_proc = 0, n_reads_last_proc = 0;

   // timers (now, at previous query)
   chrono::steady_clock::time_point time_now, time_last;
};

}

#endif
