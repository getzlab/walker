#ifndef ARGPARSE_GUARD
#define ARGPARSE_GUARD

#include <string>

namespace walker {
using namespace std;

/** A barebones set of command line arguments
 *
 */
typedef struct basic_arg {
   string bam_in;
   string output_file;
   string input_file;
   string ref_fa;
} basic_arg_t;

/** Parse basic arguments
 *  @parameter Pointer to argument struct
 *  @return Returns true if arguments parse OK, false otherwise
 */
bool basic_argparse(int argc, char** argv, basic_arg_t* args);

/** Validate basic arguments
 *  @parameter Pointer to argument struct
 *  @return Returns true if arguments validate OK, false otherwise
 */
bool basic_argparse_validate(basic_arg_t* args);

}

#endif
