#ifndef ZIGALIGN_INPUT_FASTA_H
#define ZIGALIGN_INPUT_FASTA_H

#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
using namespace std;

KSEQ_INIT(gzFile, gzread)

struct FastaSeq {
    string name;
    string seq;
};

FastaSeq input_fasta_seq(const char *fn);

double cputime();
double realtime();

#endif //ZIGALIGN_INPUT_FASTA_H
