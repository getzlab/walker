#include "walker.hpp"

// lots of this code derived from rrBAM

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
   if(!reference.LoadIndex(ref_fa)) {
      fprintf(stderr, "Error in %s : couldn't load reference!\n", __func__);
      exit(-1);
   }
}

bool walker::EDz(const SeqLib::BamRecord& record) {
   int32_t nm;
   return record.GetIntTag("NM", nm) ? nm == 0 : true;
}

vector<uint32_t> walker::nonref_pos(const SeqLib::BamRecord& record) { // {{{
   vector<uint32_t> output;

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
   for(auto c_f : c) {
      switch(c_f.Type()) {
	 /* this operator consumes bases over the ref and read
	  */
	 case 'M' :
	    for(int i = 0; i < c_f.Length(); i++) {
	       if(ror[refpos] != readseq[readpos]) {
		  output.push_back(record.PositionWithSClips() + refpos);
	       }
	       readpos++;
	       refpos++;
	    }
	    break;

	 // this operator consumes bases of the ref and read, but we don't assess
	 // their (mis)match status.
	 case 'S' :
	    readpos += c_f.Length();
	    refpos += c_f.Length();

	    break;

	 // these operators consume bases over the ref but not the read
	 case 'D' :
	 case 'N' : 
	    output.push_back(record.PositionWithSClips() + refpos);
	    refpos += c_f.Length();

	    break;

	 // this operator consumes bases over the read but not the ref
	 case 'I' :
	    output.push_back(record.PositionWithSClips() + refpos);
	    readpos += c_f.Length();

	    break;
      }
   }

   return output;
} // }}}

}
