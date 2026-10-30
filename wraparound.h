#ifndef ZIGALIGN_WRAPAROUND_H
#define ZIGALIGN_WRAPAROUND_H

#include "zigalign.h"

const int BT_VERTICAL = 1;
const int BT_HORIZONTAL = 2;
const int BT_DIAGONAL = 3;
const char UNIT_SI_CHAR = '$'; // Mark for repeat units in SI sequence

int global_alignment(const ScoreParameters &sp, int n, const char *a, int m, const char *b);

struct AlnSta {
    int mat, mis, ins, del;
    string ext_a, ext_b;
    AlnSta() {
        mat = mis = ins = del = 0;
    }
    void operator += (const AlnSta &b) {
        this->mat += b.mat;
        this->mis += b.mis;
        this->ins += b.ins;
        this->del += b.del;
    }
};

AlnSta global_cigar(const ScoreParameters &sp, int n, const char *a, int m, const char *b, int bw);

struct GlobalResult {
    int score, mat, mis, gap;
    GlobalResult() {
        score = mat = mis = gap = 0;
    }
};

GlobalResult global_cigar2(const ScoreParameters &sp, int bw, int n, const char *a, int m, const char *b);

vector<RepInterval> wraparound_dp(int pat_len, const char *pat_seq, int que_len, const char *que_seq);

struct SgResult {
    int beg, end;
    int score, mis, gap;
    SgResult() {
        beg = end = -1;
        score = -INF;
        mis = gap = 0;
    }
};

SgResult semi_global(const ScoreParameters &sc, int n, const char *a, int m, const char *b);

#endif //ZIGALIGN_WRAPAROUND_H
