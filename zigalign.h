#ifndef ZIGALIGN_ZIGALIGN_H
#define ZIGALIGN_ZIGALIGN_H

#include <vector>
#include <cstdint>
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include "input_fasta.h"
using namespace std;

const int INF = 100000000;
const int SA_MAX_LEN = 50000; // Self-alignment max length
const int SA_PART_LEN = 20000;
const double MIN_MATCH_RATIO = 0.6;
const double MAX_TR_DIV = 0.4;
const bool DEBUG = false;

struct ScoreParameters {
    int mat, mis, gap_o, gap_e;
    ScoreParameters(int ma, int mi, int go, int ge): mat(ma), mis(mi), gap_o(go), gap_e(ge) {}
};

// Zigalign parameters
struct ZigOptions {
    // SI model
    ScoreParameters si{1, -4, -6, -1};
    // DSI model
    int min_unit_size = 50;
    int open_tr_pen = -12;
    int close_tr_pen = -36;
    ScoreParameters self_aln{1, -1, -4, -1};
    ScoreParameters unit_aln{1, -120, -180, -30};
    int del_unit_pen = -20;
    int min_reps_cnt = 20;
    int max_units_del = 40000;

    int n_threads = 8;
    char cigar_neq = false; // True if cigar includes =X instead of M
    const char *log_prefix = nullptr; // Intermediate results
    string dev_label = "0000"; // Run a stage or input from local files

    // ----------------------------------------- //
    // Scoring parameters for repeats in pairwise alignment
    int mat_score;
    int mis_pen;
    int gap_o;
    int gap_e;
    // Scoring matrix for repeats in self-alignment
    int sa_mat_score;
    int sa_mis_pen;
    int sa_gap_o;
    int sa_gap_e;
    // Penalty for deleting duplications
    int del_dup_o;
    int del_dup_e;
    // K-band size
    int band_width;
    ZigOptions() {
        // Scoring parameters should be adjusted by duplication and variation rate
        mat_score = 1;
        mis_pen = -120; // Increased penalty for small variants, i.e., duplication indels are preferred.
        gap_o = -180;
        gap_e = -30;
        min_unit_size = 100;
        sa_mat_score = 2;
        sa_mis_pen = -3;
        sa_gap_o = -3;
        sa_gap_e = -1;
        open_tr_pen = -2;
        close_tr_pen = -6;
        del_unit_pen = -20;
        del_dup_o = -6;
        del_dup_e = -2;
        band_width = 50;
        n_threads = 8;
        log_prefix = nullptr;
    }
};

struct Dp1Cell {
    int D, E, F, H; // D: score for duplication events
    int beg, end; // The range of duplication or new copy
    Dp1Cell() {
        D = E = F = H = -INF;
        beg = end = -1;
    }
};

// Compact structure for backtrace
struct Bt1Cell {
    uint8_t event;
    uint8_t pi;
    uint16_t pj;
    Bt1Cell() {
        event = -1;
        pi = -1;
        pj = -1;
    }
};

struct RepInterval {
    int beg, end; // [beg, end)
    int mat, mis, gap; // #matches, mismatches and gaps

    RepInterval() {
        beg = end = -1;
        mat = mis = gap = 0;
    }

    bool operator < (const RepInterval &b) const {
        return (this->beg != b.beg) ? this->beg < b.beg : this->end < b.end;
    }
};

vector<RepInterval> self_alignment(const ZigOptions &opt, int n, const char *seq, vector<vector<Bt1Cell>> &bt);

void align_etr(const ZigOptions &opt, const char *fn_t, const char *fn_q);

#endif //ZIGALIGN_ZIGALIGN_H
