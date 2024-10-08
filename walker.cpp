#include "walker.hpp"

#include <assert.h>

// lots of this code derived from rrBAM

#define CLAMPH(x, high) (((x) > (high)) ? (high) : (x))

#define PACK_CIG(type, length) (((type) << 14) | (CLAMPH((length), 0x3FFF) & 0x3FFF))

using namespace std;

namespace walker {

walker::walker(const string& bam_in, const string& ref_fa) {
   // load BAM
   if(!this->reader.Open(bam_in)) {
      fprintf(stderr, "Error in %s : couldn't open BAM file!\n", __func__);
      exit(-1);
   }

   // load its header
   header = reader.Header();

   // load reference
   if(ref_fa != "") {
      if(!reference.LoadIndex(ref_fa)) {
         fprintf(stderr, "Error in %s : couldn't load reference!\n", __func__);
         exit(-1);
      }
   }
}

walker::~walker() {
   if(outfile != NULL) fclose(outfile);
}

bool walker::EDz(const SeqLib::BamRecord& record) {
   int32_t nm;
   return record.GetIntTag("NM", nm) ? nm == 0 : true;
}

vector<nonref_pos_t> walker::nonref_pos(const SeqLib::BamRecord& record) { // {{{
   /* lower 32 bits contain absolute reference position of mismatch
    * high 32 bits contain: 
    * - lower 16: position along the 
    * - upper 16: 2 bits for CIGAR op (as defined in sam.h, with caveat for S),
    *    14 bits for op. length. Note that CIGAR op "S" (soft clip) is 4 in sam.h,
    *    but it is 3 here (to fit in 2 bits). Since we overload 3 with "N" (RNA-seq skip),
    *    these ops must be distinguished based on their position in the read --
    *    soft clips can only occur at the read beginnings/ends, whereas "N" can never
    *    occur there.
    */
   vector<nonref_pos_t> output;

   // query reference over this read

   // for contigs with no N padding (e.g. chrM), reads can overhang the start of
   // the reference.
   // let's just ignore such reads.
   if(record.PositionWithSClips() < 0) return output;

   string ror = reference.QueryRegion(header.IDtoName(record.ChrID()), 
                                      record.PositionWithSClips(),
                                      record.PositionEndWithSClips());

   string readseq = record.Sequence();
   SeqLib::Cigar c = record.GetCigar();

   int readpos = 0;
   int refpos = 0;
   int jointpos = 0; // position over joint read/ref alignment, e.g.
   /*
      0    4    8
      ACGTACGTA--CGT REF
      ACGCAC--ACCCGT READ
         ^  ^^ ^^
         3  56 78
    */
   for(auto c_f : c) {
      switch(c_f.Type()) {
         /* this operator consumes bases over the ref and read
          */
         case 'M' :
            /*   IDEA: since we expect the set of mismatches to be sparse,
                       rather than comparing each byte at a time, can we 
                       cast pointer to int64, bitwise xor read and reference,
                       and compare 8 bytes at a time?
                       
                       could potentially be even faster if we do a binary search
                       for nonzero xor'd bits through the int64, rather than 
                       always iterating over all 8 bytes
             */
           for(int i = 0; i < c_f.Length(); i++) {
              if(pack_2bit[ror[refpos]] != pack_2bit[readseq[readpos]]) {
                 output.push_back({
                   .refpos = record.PositionWithSClips() + refpos,
                   .readpos = CLAMPH(readpos, 0xFFFF) & 0xFFFF,
                   .jointpos = CLAMPH(jointpos, 0xFFFF) & 0xFFFF,
                   .cig = PACK_CIG(0, 1)
                 });
               }
               readpos++;
               refpos++;
               jointpos++;
            }
            break;

         // this operator consumes bases of the ref and read, but we don't assess
         // their (mis)match status; we just add them to the vector.
         // note that we denote "S" as 3, which deviates from definition in sam.h
         // (see above)
         case 'S' :
            output.push_back({
              .refpos = record.PositionWithSClips() + refpos,
              .readpos = CLAMPH(readpos, 0xFFFF) & 0xFFFF,
              .jointpos = CLAMPH(jointpos, 0xFFFF) & 0xFFFF,
              .cig = PACK_CIG(3, c_f.Length())
            });
            readpos += c_f.Length();
            refpos += c_f.Length();
            jointpos += c_f.Length();

            break;

         // these operators consume bases over the ref but not the read
         case 'D' :
         case 'N' : 
            output.push_back({
              .refpos = record.PositionWithSClips() + refpos,
              .readpos = CLAMPH(readpos, 0xFFFF) & 0xFFFF,
              .jointpos = CLAMPH(jointpos, 0xFFFF) & 0xFFFF,
              .cig = PACK_CIG(c_f.RawType(), c_f.Length())
            });
            refpos += c_f.Length();
            jointpos += c_f.Length();

            break;

         // this operator consumes bases over the read but not the ref
         case 'I' :
            output.push_back({
              .refpos = record.PositionWithSClips() + refpos,
              .readpos = CLAMPH(readpos, 0xFFFF) & 0xFFFF,
              .jointpos = CLAMPH(jointpos, 0xFFFF) & 0xFFFF,
              .cig = PACK_CIG(c_f.RawType(), c_f.Length())
            });
            readpos += c_f.Length();
            jointpos += c_f.Length();

            break;
      }
   }

   return output;
} // }}}

//
// default iterators {{{

void walker::walk() {
   while(reader.GetNextRecord(cur_read)) {
      // get stats
      if(!(n_reads % 100000)) print_status();
      n_reads++;

      // apply default filters
      if(filter_read(cur_read)) continue;

      // apply user-specified action to cur_read
      if(!walk_apply(cur_read)) break;
   }
}

void walker::walk(const SeqLib::GenomicRegion& region) {
   if(!reader.SetRegion(region)) {
      fprintf(stderr, "Error setting region!\n");
      exit(1);
   }

   walk();
}

void walker::walk(const SeqLib::GenomicRegionCollection<>& region_collection) {
   if(!reader.SetMultipleRegions(region_collection)) {
      fprintf(stderr, "Error setting multiple regions!\n");
      exit(1);
   }
   walk();
}

// }}}

void walker::increment_pos(uint16_t& curchr, uint32_t& curpos) {
   // TODO: ensure that curchr isn't out of bounds
   if(curpos >= header.GetSequenceLength(curchr) - 1) {
      curchr++;
      curpos = 0;
   } else curpos++;
}

bool walker::filter_read(const SeqLib::BamRecord& record) {
   // skip if this read is unmapped
   if(!record.MappedFlag()) return true;

   // skip if this read is vendor fail
   if(record.QCFailFlag()) return true;

   // skip if read is a duplicate
   if(record.DuplicateFlag()) return true;

   // skip if read is not primary
   if(record.shared_pointer().get()->core.flag & 0x900) return true;

   // skip if read is MQZ
   if(record.MapQuality() == 0) return true;

   return false;
}

void walker::print_status() {
   if(n_reads == 0) {
      time_last = chrono::steady_clock::now();
      return;
   }

   // count reads per second
   time_now = chrono::steady_clock::now();
   double RPS = ((double) (n_reads - n_reads_last)) /
                chrono::duration_cast<chrono::duration<double>>(time_now - time_last).count();
   n_reads_last = n_reads;
   double RPS_proc = ((double) (n_reads_proc - n_reads_last_proc)) /
                chrono::duration_cast<chrono::duration<double>>(time_now - time_last).count();
   n_reads_last_proc = n_reads_proc;
   time_last = chrono::steady_clock::now();

   // print status
   if(n_reads_proc > 0) fprintf(stderr,
     "%s:%d (%0.2f, %0.2f r/s [tot., proc.])\n",
     header.IDtoName(cur_read.ChrID()).c_str(),
     cur_read.Position(),
     RPS,
     RPS_proc
   );
   else fprintf(stderr,
     "%s:%d (%0.2f r/s)\n",
     header.IDtoName(cur_read.ChrID()).c_str(),
     cur_read.Position(),
     RPS
   );

   return;
}

bool walker::set_output_file(const string& outfile) {
   if(outfile == "-") this->outfile = stdout;
   else {
      this->outfile = fopen(outfile.c_str(), "w");
      if(!this->outfile) {
         fprintf(stderr, "Could not open output file for writing!\n");
         return 0;
      }
      this->outfile_name = outfile;
   }

   return 1;
}

bool walker::close_output_file() {
   assert(fileno(this->outfile) > 2); // should not be closing standard I/O streams
   if(fclose(this->outfile)) return 0;
   this->outfile = NULL;
   return 1;
}

}
