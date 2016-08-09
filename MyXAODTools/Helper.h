#ifndef __MYXAODTOOLS_HELPER_H__
#define __MYXAODTOOLS_HELPER_H__

#include <TChain.h>

using namespace std;
namespace MyXAODTools{
    namespace Helper{
        TChain* loader(const char* inFile_name, const char* chain_name = "physics");
    }
}
#endif
