#include "argparse.hpp"

#include <iostream>

#include <sys/stat.h>
#include <unistd.h>

namespace walker {

using namespace std;

bool basic_argparse(int argc, char** argv, basic_arg_t* args) {
   // placeholder struct for arguments
   basic_arg_t args_ph = {
     .bam_in = args->bam_in != "" ? args->bam_in : "",
     .output_file = args->output_file != "" ? args->output_file : "-", // default to stdout
     .input_file = args->input_file != "" ? args->input_file : "",
     .ref_fa = args->ref_fa != "" ? args->ref_fa : ""
   };

   char arg;
   opterr = 0;
   while((arg = getopt(argc, argv, "-b:o:i:r:")) != -1) {
      switch(arg) {
         case 'b' : // path to BAM
            args_ph.bam_in = string(optarg);
            break;
         case 'o' : // output filename
            args_ph.output_file = string(optarg);
            break;
         case 'i' : // input sites list
            args_ph.input_file = string(optarg);
            break;
         case 'r' : // path to reference
            args_ph.ref_fa = string(optarg);
            break;
         case '?' :
            string missing;
            switch(optopt) {
               case 'b' :
                  missing = "Path to BAM";
                  break;
               case 'o' :
                  missing = "Output filename";
                  break;
               case 'i' :
                  missing = "Input filename";
                  break;
               case 'r' :
                  missing = "Path to reference";
                  break;

             cerr << missing + " cannot be blank!\n";
             return false;
          }
      }
   }

   // we wait until getopt finishes successfully before populating the external argument struct
   *args = args_ph;

   return true;
}

bool basic_argparse_validate(basic_arg_t* args) {
   // will check if BAM is openable later.

   // if BAM is a GCS bucket path, make sure the OAUTH token is set in the environment
   /*if(bam_in.compare(0, 3, "gs:") == 0 && getenv("GCS_OAUTH_TOKEN") == NULL) {
      cerr << "Error: GCS OAUTH token not set!\n";
      return 1;
   }*/

   // check that reference FASTA is present and a regular file
   struct stat sb;
   if(stat(args->ref_fa.c_str(), &sb) == -1 || !S_ISREG(sb.st_mode)) {
      cerr << "Error: invalid reference FASTA!\n";
      return false;
   }

   return true;
}

}
