#include <cmath>
#include <memory>
#include <sstream>

#include "zigalign.h"
#include "wraparound.h"
#include "union_set.h"

vector<RepInterval> self_alignment(const ZigOptions &opt, int n, const char *seq, vector<vector<Bt1Cell>> &bt)
{
	if (n > SA_MAX_LEN) {
		fprintf(stderr, "Sequence length %d exceeds the limit of self alignment %d\n", n, SA_MAX_LEN);
		abort();
	}

	const int MIN_UNIT = opt.min_unit_size;
	const int MAT_SCORE = 0;
	const int OPEN_TR = opt.open_tr_pen;
	const int CLOSE_TR = opt.close_tr_pen;
	const int SA_MAT_SCORE = opt.self_aln.mat;
	const int SA_MIS_PEN = opt.self_aln.mis;
	const int SA_GAP_O = opt.self_aln.gap_o;
	const int SA_GAP_E = opt.self_aln.gap_e;

	const int BT_NON_REP = 0;
	const int BT_START_REP = 1;
	const int BT_NEW_COPY = 2;
	const int BT_WITHIN_REP = 3;
	const int BT_END_REP = 4;

	if (MAT_SCORE >= SA_MAT_SCORE) {
		fprintf(stderr, "Warning: setting higher score for diagonal match can't detect tandem repeats\n");
	}

	if (n < MIN_UNIT) {
		vector<RepInterval> ret;
		return ret;
	}

	// Memory allocation
	vector<Dp1Cell> prev_dp(n + 1);
	vector<Dp1Cell> curr_dp(n + 1);
 	bt.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		bt[i].resize(i + 1);
	}

	// Initialization
	prev_dp[0].H = 0;
	prev_dp[0].beg = 1;
	for (int i = 1; i <= MIN_UNIT; i++) {
		bt[i][i].event = BT_NON_REP;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
	}
	prev_dp[MIN_UNIT].H = MIN_UNIT * MAT_SCORE;
	prev_dp[MIN_UNIT].beg = 1;

	// Main loop
	for (int i = MIN_UNIT + 1; i <= n; i++) {
		for (int j = 0; j <= i; j++) curr_dp[j] = Dp1Cell();
		curr_dp[i].H = prev_dp[i-1].H + MAT_SCORE;
		bt[i][i].event = BT_NON_REP;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
		curr_dp[i].beg = prev_dp[i-1].beg;

		int max_value = (i-1) * MAT_SCORE; // From non-repetitive region (NOTE: decreased score)
		int t_beg = prev_dp[i-1].beg; // To prevent illegal path
		int t_end = i - 1;
		uint8_t event = BT_START_REP;
		for (int j = 1; j < i-1; j++) {
			// Start new copy at the end of duplications
			if (prev_dp[j].end == j and prev_dp[j].H > max_value) {
				max_value = prev_dp[j].H;
				t_beg = prev_dp[j].beg;
				t_end = j;
				event = BT_NEW_COPY;
			}
		}

		// FIXME: maximum value should be chosen from overlapping template intervals
		// D transfer
		for (int j = t_beg; j <= t_end - MIN_UNIT + 1; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? SA_MAT_SCORE : SA_MIS_PEN;
			int pen = (event == BT_START_REP) ? OPEN_TR : 0; // No penalty for new copy
			curr_dp[j].D = max_value + tmp + pen;
			curr_dp[j].beg = j;
			curr_dp[j].end = t_end;
			bt[i][j].event = event;
			bt[i][j].pi = 1;
			bt[i][j].pj = t_end;
		}
		for (int j = t_end - MIN_UNIT + 2; j <= i - MIN_UNIT; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? SA_MAT_SCORE : SA_MIS_PEN;
			curr_dp[j].D = (i - 1) * MAT_SCORE + tmp + OPEN_TR;
			curr_dp[j].beg = j;
			curr_dp[j].end = i - 1;
			bt[i][j].event = BT_START_REP;
			bt[i][j].pi = 1;
			bt[i][j].pj = i - 1;
		}

		// Sub matrix of repetition
		for (int j = 1; j < i; j++) {
			int v_score = -INF, h_score = -INF, d_score = -INF;
			if (j >= prev_dp[j].beg and j <= prev_dp[j].end) {
				v_score = max(max(prev_dp[j].D, prev_dp[j].H) + SA_GAP_O, prev_dp[j].E) + SA_GAP_E;
				curr_dp[j].E = v_score;
			}
			if (j >= curr_dp[j-1].beg and j <= curr_dp[j-1].end) {
				h_score = max(max(curr_dp[j-1].D, curr_dp[j-1].H) + SA_GAP_O, curr_dp[j-1].F) + SA_GAP_E;
				curr_dp[j].F = h_score;
			}
			if (j >= prev_dp[j-1].beg and j <= prev_dp[j-1].end) {
				int tmp = (seq[i-1] == seq[j-1] ? SA_MAT_SCORE : SA_MIS_PEN);
				d_score = max(prev_dp[j-1].D, prev_dp[j-1].H) + tmp;
			}

			// Here, >= prefers continue an existing copy, so it can generate longer tandem repeats
			if (v_score >= curr_dp[j].D and v_score > curr_dp[j].H) {
				curr_dp[j].H = v_score;
				curr_dp[j].beg = prev_dp[j].beg;
				curr_dp[j].end = prev_dp[j].end;
				bt[i][j].event = BT_WITHIN_REP;
				bt[i][j].pi = 1;
				bt[i][j].pj = j;
			}
			if (h_score >= curr_dp[j].D and h_score > curr_dp[j].H) {
				curr_dp[j].H = h_score;
				curr_dp[j].beg = curr_dp[j-1].beg;
				curr_dp[j].end = curr_dp[j-1].end;
				bt[i][j].event = BT_WITHIN_REP;
				bt[i][j].pi = 0;
				bt[i][j].pj = j-1;
			}
			if (d_score >= curr_dp[j].D and d_score > curr_dp[j].H) {
				curr_dp[j].H = d_score;
				curr_dp[j].beg = prev_dp[j-1].beg;
				curr_dp[j].end = prev_dp[j-1].end;
				bt[i][j].event = BT_WITHIN_REP;
				bt[i][j].pi = 1;
				bt[i][j].pj = j-1;
			}
		}

		// B transfer
		max_value = -INF;
		t_beg = -1;
		t_end = -1;
		for (int j = 1; j < i; j++) {
			// Only return to diagonal if sub-matrix reaches the lower-right corner
			if (curr_dp[j].H > max_value and curr_dp[j].end == j) {
				max_value = curr_dp[j].H;
				t_beg = curr_dp[j].beg;
				t_end = curr_dp[j].end;
			}
		}
		if (max_value + CLOSE_TR > curr_dp[i].H) {
			curr_dp[i].H = max_value + CLOSE_TR;
			curr_dp[i].beg = t_beg; // This variable is reused to prevent illegal path
			bt[i][i].event = BT_END_REP;
			bt[i][i].pi = 0;
			bt[i][i].pj = t_end;
		}

		swap(prev_dp, curr_dp);
	}

	// Backtrace
	int ti = n, tj = n;
	vector<RepInterval> reps;
	while (ti > 0 and tj > 0) {
		Bt1Cell t = bt[ti][tj];
		if (t.event == BT_END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			while (t.event != BT_START_REP) {
				ti = t.pi == 1 ?ti-1 :ti;
				tj = t.pj;
				t = bt[ti][tj]; // Lower-right corner of the sub-matrix
				RepInterval u;
				u.end = ti + 1;
				while (t.event != BT_NEW_COPY and t.event != BT_START_REP) {
					if (t.pi == 1 and tj == t.pj + 1) {
						u.mis += (seq[ti - 1] != seq[tj - 1]);
						u.mat += (seq[ti - 1] == seq[tj - 1]);
					} else {
						u.gap++;
					}
					ti = t.pi == 1 ?ti-1 :ti;
					tj = t.pj;
					t = bt[ti][tj];
				}
				// Pointer is now at the upper-left corner
				u.mis += (seq[ti - 1] != seq[tj - 1]);
				u.mat += (seq[ti - 1] == seq[tj - 1]);
				u.beg = ti;
				reps.push_back(u);
			}
			// Template for all copy units above
			RepInterval tem;
			tem.beg = tj;
			tem.end = ti;
			tem.mis = tem.gap = -1;
			reps.push_back(tem);
		}
		ti = t.pi == 1 ?ti-1 :ti;
		tj = t.pj;
	}

	if (reps.empty()) return reps;

	// Reset to 0-based index
	for (RepInterval &r: reps) {
		r.beg--;
		r.end--;
	}

	// de-overlap repeat units (tandem repeats should be neatly stacked)
	sort(reps.begin(), reps.end());
	int new_size = 0;
	for (int i = 0; i < reps.size(); ) {
		const RepInterval &x = reps[i];
		bool kept = true;
		int j = i + 1;
		for (; j < reps.size(); j++) {
			const RepInterval &y = reps[j];
			if (y.beg >= x.end) break;
			// Keep the longer one
			if (y.end - y.beg >= x.end - x.beg) {
				kept = false;
				break;
			}
		}
		if (kept) reps[new_size++] = x;
		i = j;
	}
	reps.resize(new_size);
	return reps;
}

struct TandemRepeats {
	int pid; // Pattern ID
	int base_idx; // Base index
	double ave_len, ave_mis, ave_gap; // Meta information
	vector<RepInterval> repeats;
	TandemRepeats():pid(-1), base_idx(0), ave_len(0), ave_mis(0), ave_gap(0) {}
};

vector<TandemRepeats> group_tandem_repeats(const vector<RepInterval> &reps)
{
	vector<TandemRepeats> ret;
	if (reps.empty()) return ret;

	TandemRepeats tr;
	tr.ave_len = reps[0].end - reps[0].beg;
	tr.ave_mis = tr.ave_gap = 0;
	tr.repeats.push_back(reps[0]);
	for (int i = 1; i < reps.size(); i++) {
		const RepInterval &r = reps[i];
		if (r.mat == -1) { // New group
			tr.ave_len /= tr.repeats.size();
			tr.ave_mis /= (tr.repeats.size() - 1);
			tr.ave_gap /= (tr.repeats.size() - 1);
			ret.push_back(tr);
			tr.ave_len = r.end - r.beg;
			tr.ave_mis = tr.ave_gap = 0;
			tr.repeats.push_back(r);
		} else {
			tr.ave_len += r.end - r.beg;
			tr.ave_gap += r.gap;
			tr.ave_mis += r.mis;
			tr.repeats.push_back(r);
		}
	}
	tr.ave_len /= tr.repeats.size();
	tr.ave_mis /= (tr.repeats.size() - 1);
	tr.ave_gap /= (tr.repeats.size() - 1);
	ret.push_back(tr);
	return ret;
}

// ---- for debug ----
void output_repeats(const string &fn, const vector<TandemRepeats> &lr)
{
	ofstream out(fn);
	assert(out.is_open());
	size_t n = lr.size();
	out.write((char*)&n, sizeof(size_t));
	for (const TandemRepeats &tr: lr) {
		out.write((char*) &tr.pid, sizeof(int));
		out.write((char*) &tr.base_idx, sizeof(int));
		out.write((char*) &tr.ave_len, sizeof(double));
		out.write((char*) &tr.ave_mis, sizeof(double));
		out.write((char*) &tr.ave_gap, sizeof(double));
		n = tr.repeats.size();
		out.write((char*) &n, sizeof(size_t));
		for (const RepInterval &r: tr.repeats) {
			out.write((char*) &r.beg, sizeof(int));
			out.write((char*) &r.end, sizeof(int));
			out.write((char*) &r.mat, sizeof(int));
			out.write((char*) &r.mis, sizeof(int));
			out.write((char*) &r.gap, sizeof(int));
		}
	}
	out.close();
}

vector<TandemRepeats> input_repeats(const string &fn)
{
	vector<TandemRepeats> ret;
	ifstream in(fn);
	assert(in.is_open());
	size_t n;
	in.read((char*)&n, sizeof(size_t));
	ret.resize(n);
	for (int i = 0; i < n; i++) {
		TandemRepeats &tr = ret[i];
		in.read((char*)&tr.pid, sizeof(int));
		in.read((char*)&tr.base_idx, sizeof(int));
		in.read((char*)&tr.ave_len, sizeof(double));
		in.read((char*)&tr.ave_mis, sizeof(double));
		in.read((char*)&tr.ave_gap, sizeof(double));

		size_t m;
		in.read((char*)&m, sizeof(size_t));
		tr.repeats.resize(m);
		for (int j = 0; j < m; j++) {
			RepInterval &r = tr.repeats[j];
			in.read((char*) &r.beg, sizeof(int));
			in.read((char*) &r.end, sizeof(int));
			in.read((char*) &r.mat, sizeof(int));
			in.read((char*) &r.mis, sizeof(int));
			in.read((char*) &r.gap, sizeof(int));
		}
	}
	in.close();

	for (int i = 0; i < ret.size(); i++) {
		TandemRepeats &tr = ret[i];
		fprintf(stderr, "TR %d: size=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			i + 1, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	fprintf(stderr, "\n");
	return ret;
}

// Test if group B can be merged into group A
bool merge_tr_test(const TandemRepeats &tr_a, const char *seq_a, const TandemRepeats &tr_b, const char *seq_b)
{
	// Average length is more reliable
	double len_diff = abs(tr_a.ave_len - tr_b.ave_len) / tr_a.ave_len;
	if (len_diff > MAX_TR_DIV) return false;

	// Use average divergence to reduce the disturbance of outliers
	const int MAX_COMPARE = 20;
	int step_a = (tr_a.repeats.size() + MAX_COMPARE - 1) / MAX_COMPARE;
	int step_b = (tr_b.repeats.size() + MAX_COMPARE - 1) / MAX_COMPARE;
	int n_comp = 0;
	double ave_div = 0;
	for (int i = 0; i < tr_b.repeats.size(); i += step_b) {
		const RepInterval &b = tr_b.repeats[i];
		int len_b = b.end - b.beg;
		string sb;
		sb.resize(len_b);
		memcpy((char*) sb.data(), seq_b + b.beg, len_b);
		sb += sb;
		for (int j = 0; j < tr_a.repeats.size(); j += step_a) {
			const RepInterval &a = tr_a.repeats[j];
			int len_a = a.end - a.beg;
			// Use semi-global alignment to handle repeat rotations
			SgResult sg = semi_global(len_a, seq_a + a.beg, sb.length(), sb.data());
			n_comp++;
			ave_div += 1.0 * (sg.mis + sg.gap) / len_a;
		}
	}
	ave_div /= n_comp;
	return ave_div < MAX_TR_DIV;
}

vector<TandemRepeats> collect_long_repeats(const ZigOptions &opt, int n, const char *seq)
{
	int n_threads = opt.n_threads;
	vector<vector<Bt1Cell>> bt_bin[n_threads];
	vector<RepInterval> reps_bin[n_threads];
	vector<TandemRepeats> tr_bin[n_threads];
	vector<TandemRepeats> ret;
	int global_os = 0;
	int batch_id = 0;
	while (global_os < n) {
		double r_start = realtime(), c_start = cputime();
		int m = n - global_os;
		const char *s0 = seq + global_os;

		// Multi-threading DP
		#pragma omp parallel for
		for (int i = 0; i < n_threads; i++) {
			int beg = i * SA_PART_LEN;
			int end = min(m, beg + SA_PART_LEN);
			if (beg >= end) {
				reps_bin[i].clear();
				tr_bin[i].clear();
			} else {
				const char *s = s0 + beg;
				vector<RepInterval> curr_reps, next_reps;
				curr_reps = self_alignment(opt, end - beg, s, bt_bin[i]);
				reps_bin[i].clear();
				while (not curr_reps.empty()) {
					next_reps.clear();
					for (const RepInterval &r: curr_reps) {
						// fprintf(stderr, "Break down [%d, %d) len=%d\n", r.beg, r.end, r.end - r.beg);
						vector<RepInterval> tmp = self_alignment(opt, r.end - r.beg, s + r.beg, bt_bin[i]);
						if (tmp.empty()) {
							// Non-nested
							reps_bin[i].push_back(r);
						} else {
							// Recursively break down the nested repeat
							for (RepInterval &t: tmp) {
								t.beg += r.beg;
								t.end += r.beg;
								next_reps.push_back(t);
							}
						}
					}
					swap(curr_reps, next_reps);
				}
				sort(reps_bin[i].begin(), reps_bin[i].end());
				// Remove the last unit because it might be incomplete
				if (not reps_bin[i].empty()) {
					reps_bin[i].resize(reps_bin[i].size() - 1);
				}
				for (RepInterval &r: reps_bin[i]) {
					r.beg += global_os + beg;
					r.end += global_os + beg;
				}
				tr_bin[i] = group_tandem_repeats(reps_bin[i]);
			}
		}
		double c_time = cputime() - c_start, r_time = realtime() - r_start;
		int n_reps = 0;
		for (int i = 0; i < n_threads; i++) n_reps += reps_bin[i].size();
		fprintf(stderr, "Batch %d, offset: %d, reps: %d; CPU: %.2f, Wall: %.2f, Core: %.2f\n",
			++batch_id, global_os, n_reps, c_time, r_time, c_time / r_time);

		// Bridge adjacent batches
		for (int i = 0; i < n_threads; i++) {
			for (const TandemRepeats &g: tr_bin[i]) {
				if (ret.empty()) {
					ret.push_back(g);
				} else {
					TandemRepeats &t = ret.back();
					bool to_merge = merge_tr_test(t, seq, g, seq);
					if (to_merge) {
						t.repeats.insert(t.repeats.end(), g.repeats.begin(), g.repeats.end());
					} else {
						ret.push_back(g);
					}
				}
			}
		}
		global_os += n_threads * SA_PART_LEN;
	}

	int new_n = 0;
	for (int i = 0; i < ret.size();i ++) {
		TandemRepeats &tr = ret[i];
		// Is it correct to skip the small tandem repeats?
		if (tr.repeats.size() < opt.min_reps_cnt) continue;
		tr.pid = -1;
		tr.base_idx = 0;
		tr.ave_len = tr.ave_gap = tr.ave_mis = 0;
		int m = 0;
		for (const RepInterval &r: tr.repeats) {
			tr.ave_len += r.end - r.beg;
			if (r.mat != -1) {
				m++;
				tr.ave_gap += r.gap;
				tr.ave_mis += r.mis;
			}
		}
		tr.ave_len /= tr.repeats.size();
		if (m) tr.ave_mis /= m;
		if (m) tr.ave_gap /= m;
		ret[new_n++] = tr;
	}
	ret.resize(new_n);

	for (int i = 0; i < ret.size(); i++) {
		TandemRepeats &tr = ret[i];
		fprintf(stderr, "TR %d: size=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			i + 1, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	fprintf(stderr, "\n");

	return ret;
}

void extended_paf_format(const string &t_name, int t_len, const string &ext_t,
	const string &q_name, int q_len, const string ext_q, const vector<int> &rid_map)
{
	assert(ext_t.length() == ext_q.length());
	char strand = '+';
	int matches_n = 0, mismatches_n = 0, len = ext_t.length();
	int mapq = 60;

	string cigar;
	char op_type = '@';
	int op_cnt = 0;
	int t_pointer = 0;
	for (int i = 0; i < len; i++) {
		if (ext_t[i] == 'A' or ext_t[i] == 'C' or ext_t[i] == 'G' or ext_t[i] == 'T') {
			t_pointer++;
		}
		if (ext_t[i] == '[' or ext_t[i] == ']') {
			if (op_type != '@') {
				cigar += to_string(op_cnt);
				cigar += op_type;
			}
			if (ext_t[i] == '[') {
				int id = -1;
				for (int j = 0; j < 20; j++) {
					if (t_pointer - j >= 0 and rid_map[t_pointer - j] != -1) {
						id = rid_map[t_pointer - j];
						break;
					}
					if (t_pointer + j < t_len and rid_map[t_pointer + j] != -1) {
						id = rid_map[t_pointer + j];
						break;
					}
				}
				if (id != -1) cigar += '(' + to_string(id) + ')';
			}
			cigar += ext_t[i];
			op_type = '@';
			op_cnt = 0;
		} else {
			int curr_ot;
			if (ext_t[i] == '-') {
				curr_ot = 'I';
			} else if (ext_q[i] == '-') {
				curr_ot = 'D';
			} else {
				curr_ot = 'M'; // Match and mismatches
				if (ext_t[i] == ext_q[i]) matches_n++;
				else mismatches_n++;
			}

			if (curr_ot == op_type) op_cnt++;
			else {
				if (op_type != '@') {
					cigar += to_string(op_cnt);
					cigar += op_type;
				}
				op_type = curr_ot;
				op_cnt = 1;
			}
		}
	}
	if (op_type != '@' and op_cnt > 0) {
		cigar += to_string(op_cnt);
		cigar += op_type;
	}

	fprintf(stdout, "%s\t%d\t%d\t%d\t%c\t", q_name.c_str(), q_len, 0, q_len, strand);
	fprintf(stdout, "%s\t%d\t%d\t%d\t", t_name.c_str(), t_len, 0, t_len);
	fprintf(stdout, "%d\t%d\t%d\t", matches_n + mismatches_n, len, mapq);
	fprintf(stdout, "cg:Z:%s\n", cigar.c_str());
}

void print_tr_by_class(int n, const vector<TandemRepeats> &lr_t, const vector<TandemRepeats> &lr_q)
{
	vector<pair<bool, int>> group[n];
	for (int i = 0; i < lr_t.size(); i++) {
		int pid = lr_t[i].pid;
		group[pid].emplace_back(true, i);
	}
	for (int i = 0; i < lr_q.size(); i++) {
		int pid = lr_q[i].pid;
		group[pid].emplace_back(false, i);
	}
	for (int i = 0; i < n; i++) {
		if (group[i].empty()) continue;
		fprintf(stderr, "Class %d:\n", i + 1);
		for (const pair<bool, int> &p: group[i]) {
			if (p.first) {
				const TandemRepeats &tr = lr_t[p.second];
				fprintf(stderr, "  [T] reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
					tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
			} else {
				const TandemRepeats &tr = lr_q[p.second];
				fprintf(stderr, "  [Q] reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
					tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
			}
		}
	}
	fprintf(stderr, "\n");
}

void print_tr_by_order(const vector<TandemRepeats> &lr_t, const vector<TandemRepeats> &lr_q)
{
	for (const TandemRepeats &tr: lr_t) {
		fprintf(stderr, "  [T] pid=%d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			tr.pid, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	for (const TandemRepeats &tr: lr_q) {
		fprintf(stderr, "  [Q] pid=%d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			tr.pid, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	fprintf(stderr, "\n");
}

int update_pid(vector<TandemRepeats> &lr_t, vector<TandemRepeats> &lr_q)
{
	int n = 0;
	for (const TandemRepeats &tr: lr_t) n = max(n, tr.pid);
	for (const TandemRepeats &tr: lr_q) n = max(n, tr.pid);
	n++;

	int m = 0;
	vector<int> gid(n, -1);
	for (int i = 0; i < lr_t.size(); i++) {
		int pid = lr_t[i].pid;
		if (gid[pid] == -1) gid[pid] = m++;
		lr_t[i].pid = gid[pid];
	}
	for (int i = 0; i < lr_q.size(); i++) {
		int pid = lr_q[i].pid;
		if (gid[pid] == -1) gid[pid] = m++;
		lr_q[i].pid = gid[pid];
	}
	return m;
}

int categorize_repeats(const ZigOptions &opt, vector<TandemRepeats> &lr_t, const char *seq_t, vector<TandemRepeats> &lr_q, const char *seq_q)
{
	double t_real = realtime(), t_cpu = cputime();
	int n = 0;
	for (TandemRepeats &tr: lr_t) tr.pid = n++;
	for (TandemRepeats &tr: lr_q) tr.pid = n++;

	vector<pair<int,int>> tasks;
	for (int i = 0; i < lr_t.size(); i++) {
		for (int j = 0; j < lr_q.size(); j++) {
			tasks.emplace_back(i, j);
		}
	}

	int n_threads = opt.n_threads;
	vector<pair<int, int>> link_bin[n_threads];
	#pragma omp parallel for
	for (int i = 0; i < tasks.size(); i++) {
		const pair<int,int> &p = tasks[i];
		const TandemRepeats &t = lr_t[p.first];
		const TandemRepeats &q = lr_q[p.second];
		bool to_merge = merge_tr_test(t, seq_t, q, seq_q) and
						merge_tr_test(q, seq_q, t, seq_t);
		if (to_merge) {
			int tid = omp_get_thread_num();
			link_bin[tid].emplace_back(t.pid, q.pid);
		}
	}

	UnionFind us(n);
	for (int i = 0; i < n_threads; i++) {
		for (const pair<int,int> &p: link_bin[i]) {
			us.unite(p.first, p.second);
		}
	}
	for (TandemRepeats &tr: lr_t) tr.pid = us.find(tr.pid);
	for (TandemRepeats &tr: lr_q) tr.pid = us.find(tr.pid);

	// Fill in the empty slots of pattern ID
	int m = update_pid(lr_t, lr_q);
	fprintf(stderr, "Categorize %d TR groups to %d classes in %.2f CPU sec, %.2f real sec\n",
		n, m, cputime() - t_cpu, realtime() - t_real);
	return m;
}

// Picking the sanest pattern is important for normalization
string pickout_pattern(const TandemRepeats &tr, const char *seq)
{
	double min_len_diff = INF;
	RepInterval pat_rep;
	for (const RepInterval &r: tr.repeats) {
		double len_diff = abs((r.end - r.beg) - tr.ave_len);
		if (len_diff < min_len_diff) {
			min_len_diff = len_diff;
			pat_rep = r;
		}
	}
	string ret;
	int len = pat_rep.end - pat_rep.beg;
	ret.resize(len);
	memcpy((char*) ret.data(), seq + pat_rep.beg, len);
	return ret;
}

int normalize_repeats(const ZigOptions &opt, int n_class,
	int t_len, const char *t_seq, vector<TandemRepeats> &lr_t,
	int q_len, const char *q_seq, vector<TandemRepeats> &lr_q)
{
	vector<string> patterns(n_class);
	for (const TandemRepeats &tr: lr_t) {
		int pid = tr.pid;
		if (patterns[pid].empty()) {
			patterns[pid] = pickout_pattern(tr, t_seq);
		}
	}
	for (const TandemRepeats &tr: lr_q) {
		int pid = tr.pid;
		if (patterns[pid].empty()) {
			patterns[pid] = pickout_pattern(tr, q_seq);
		}
	}

	double c_time = cputime(), r_time = realtime();
	const int FLANK_EXT = 10;
	#pragma omp parallel for
	for (int k = 0; k < lr_t.size() + lr_q.size(); k++) {
		TandemRepeats &tr = k < lr_t.size() ?lr_t[k] :lr_q[k - lr_t.size()];
		const char *seq = k < lr_t.size() ?t_seq :q_seq;
		int len = k < lr_t.size() ?t_len :q_len;

		const string &pat = patterns[tr.pid];
		int l_range = tr.repeats.front().beg;
		int r_range = tr.repeats.back().end;
		// This is for the missed flanking repeats
		l_range = max(0, l_range - (int)round(pat.length() * FLANK_EXT));
		r_range = min(len, r_range + (int)round(pat.length() * FLANK_EXT));

		vector<RepInterval> reps = wraparound_dp(
			pat.length(), pat.data(), r_range - l_range, seq + l_range);
		for (RepInterval &r: reps) {
			r.beg += l_range;
			r.end += l_range;
		}
		int trim_l = 0, trim_r = 0;
		// Trim out poor flanking partitions
		for (const RepInterval &r: reps) {
			double div = 1.0 * (r.mis + r.gap) / pat.length();
			if (div > MAX_TR_DIV) {
				trim_l++;
			} else {
				break;
			}
		}
		// fprintf(stderr, "left trim %d poor partitions\n", trim_l);
		for (int i = reps.size() - 1; i >= 0; i--) {
			const RepInterval &r = reps[i];
			double div = 1.0 * (r.mis + r.gap) / pat.length();
			if (div > MAX_TR_DIV) {
				trim_r++;
			} else {
				break;
			}
		}
		// fprintf(stderr, "right trim %d poor partitions\n", trim_r);

		tr.ave_len = tr.ave_gap = tr.ave_mis = 0;
		tr.repeats.clear();
		for (int i = trim_l; i < reps.size() - trim_r; i++) {
			const RepInterval &r = reps[i];
			tr.repeats.push_back(r);
			tr.ave_len += r.end - r.beg;
			tr.ave_gap += r.gap;
			tr.ave_mis += r.mis;
		}
		// Repeats might be empty when no clear partitions emerge
		if (not tr.repeats.empty()) tr.ave_len /= tr.repeats.size();
		if (tr.repeats.size() > 1) tr.ave_gap /= (tr.repeats.size() - 1);
		if (tr.repeats.size() > 1) tr.ave_mis /= (tr.repeats.size() - 1);
	}

	// FIXME: it shouldn't happen
	for (int i = 0; i < lr_t.size() - 1; i++) {
		int beg = lr_t[i+1].repeats.front().beg;
		vector<RepInterval> &reps = lr_t[i].repeats;
		int m = reps.size();
		for (int j = reps.size()-1; j >= 0; j--) {
			const RepInterval &r = reps[j];
			if (r.end > beg) {
				m--;
			} else {
				break;
			}
		}
		reps.resize(m);
	}
	for (int i = 0; i < lr_q.size() - 1; i++) {
		int beg = lr_q[i+1].repeats.front().beg;
		vector<RepInterval> &reps = lr_q[i].repeats;
		int m = reps.size();
		for (int j = reps.size()-1; j >= 0; j--) {
			const RepInterval &r = reps[j];
			if (r.end > beg) {
				m--;
			} else {
				break;
			}
		}
		reps.resize(m);
	}

	int poor_tr = 0;
	const int MIN_REP_CNT = opt.min_reps_cnt;
	int old_n = 0;
	for (int i = 0; i < lr_t.size(); i++) {
		if (lr_t[i].repeats.size() < MIN_REP_CNT) {
			poor_tr++;
		} else {
			lr_t[old_n++] = lr_t[i];
		}
	}
	lr_t.resize(old_n);
	old_n = 0;
	for (int i = 0; i < lr_q.size(); i++) {
		if (lr_q[i].repeats.size() < MIN_REP_CNT) {
			poor_tr++;
		} else {
			lr_q[old_n++] = lr_q[i];
		}
	}
	lr_q.resize(old_n);
	n_class = update_pid(lr_t, lr_q);

	fprintf(stderr, "Partition tandem repeats in %.2f CPU sec, %.2f real sec\n",
		cputime()-c_time, realtime()-r_time);
	return n_class;
}

void collect_nested_repeats(const ZigOptions &opt, const TandemRepeats &tr, const char *seq)
{
	int n = tr.repeats.size();
	vector<int> link_t(n, -1);
	const int BACKWARD_CMP = 100;
	const double MAX_DIV = 0.1;
	#pragma omp parallel for
	for (int i = 0; i < tr.repeats.size(); i++) {
		const RepInterval &a = tr.repeats[i];
		int len_a = a.end - a.beg;
		int max_score = 0, max_id = -1;
		for (int j = max(0, i-BACKWARD_CMP); j < i; j++) {
			const RepInterval &b = tr.repeats[j];
			int len_b = b.end - b.beg;
			int bw = round(len_b * MAX_DIV);
			GlobalResult gr = global_cigar2(opt.si, bw, len_a, seq + a.beg, len_b, seq + b.beg);
			double div = 1.0 * (gr.mis + gr.gap) / min(len_a, len_b);
			if (div < MAX_DIV and gr.score > max_score) {
				max_score = gr.score;
				max_id = j;
			}
		}
		if (max_id != -1) {
			link_t[i] = max_id;
		}
	}
	vector<int> sid(n);
	for (int i = 0; i < n; i++) sid[i] = i;
	for (int i = 0; i < n; i++) {
		if (link_t[i] != -1) {
			us_union(sid, i, link_t[i]);
		}
	}
	vector<int> cnt(n, 0);
	for (int i = 0; i < n; i++) {
		int x = us_find(sid, i);
		cnt[x]++;
	}
	const int HIGH_FREQ = 100;
	const double MIN_RATIO = 0.8;
	int m = 0;
	for (int i = 0; i < n; i++) {
		int x = us_find(sid, i);
		if (x > HIGH_FREQ) {
			m++;
		}
	}
	double ratio = 1.0 * m / n;
	if (ratio > MIN_RATIO) {
		fprintf(stderr, "Observe nested structure in TR size=%ld, range=[%d,%d)\n",
			tr.repeats.size(), tr.repeats.front().beg, tr.repeats.back().end);
		for (int i = 0; i < n; i++) {
			fprintf(stderr, "%d -> %d\n", i, us_find(sid, i));
		}
		fprintf(stderr, "\n");
	}
}

vector<vector<vector<int>>> build_score_matrix(const ZigOptions &opt, int n_class,
	vector<TandemRepeats> &lr_t, int t_len, const char *t_seq,
	vector<TandemRepeats> &lr_q, int q_len, const char *q_seq)
{
	int n = 0, m = 0;
	int t_sum = 0, q_sum = 0;
	vector<vector<RepInterval>> reps_t(n_class);
	vector<vector<RepInterval>> reps_q(n_class);
	for (TandemRepeats &lr: lr_t) {
		lr.base_idx = n;
		n += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			reps_t[lr.pid].push_back(r);
			t_sum += r.end - r.beg;
		}
	}
	for (TandemRepeats &lr: lr_q) {
		lr.base_idx = m;
		m += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			reps_q[lr.pid].push_back(r);
			q_sum += r.end - r.beg;
		}
	}
	fprintf(stderr, "T: %d repeats, fraction=%.2f %%\n", n, 100.0 * t_sum / t_len);
	fprintf(stderr, "Q: %d repeats, fraction=%.2f %%\n", m, 100.0 * q_sum / q_len);

	// Build score matrix of repeat units
	double t_cpu = cputime(), t_real = realtime();
	const int MAX_UNITS_DEL = opt.max_units_del;
	vector<vector<vector<int>>> matrix(n_class);
	for (int pid = 0; pid < n_class; pid++) {
		n = reps_t[pid].size();
		m = reps_q[pid].size();
		matrix[pid].resize(n);
		for (int i = 0; i < n; i++) matrix[pid][i].resize(m, -INF);

		#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			int beg = max(i - MAX_UNITS_DEL, 0);
			int end = min(i + MAX_UNITS_DEL, m);
			const RepInterval &t = reps_t[pid][i];
			for (int j = beg; j < end; j++) {
				const RepInterval &q = reps_q[pid][j];
				matrix[pid][i][j] = global_alignment(
					opt.unit_aln, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
			}
		}
	}
	fprintf(stderr, "Build score matrix in %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);
	return matrix;
}

void output_matrix(const string &fn, const vector<vector<vector<int>>> &matrix)
{
	ofstream out(fn, ios::binary);
	assert(out.is_open());
	size_t n_class = matrix.size();
	out.write((char*)&n_class, sizeof(int64_t));
	for (int pid = 0; pid < n_class; pid++) {
		size_t n = matrix[pid].size();
		size_t m = n == 0 ?0: matrix[pid][0].size();
		out.write((char*)&n, sizeof(int64_t));
		out.write((char*)&m, sizeof(int64_t));
		if (n == 0 or m == 0) continue;
		fprintf(stderr, "pid=%d, n=%ld, m=%ld\n", pid, n, m);
		for (int i = 0; i < n; i++) {
			assert(matrix[pid][i].size() == m);
			out.write((char*)matrix[pid][i].data(), sizeof(int) * m);
		}
	}
	out.close();
}

vector<vector<vector<int>>> input_matrix(const string &fn)
{
	ifstream in(fn, ios::binary);
	assert(in.is_open());
	vector<vector<vector<int>>> matrix;
	size_t n_class;
	in.read((char*)&n_class, sizeof(int64_t));
	matrix.resize(n_class);
	fprintf(stderr, "%ld\n", n_class);
	for (int pid = 0; pid < n_class; pid++) {
		size_t n, m;
		in.read((char*)&n, sizeof(int64_t));
		in.read((char*)&m, sizeof(int64_t));
		if (n == 0 or m == 0) continue;
		fprintf(stderr, "pid=%d, n=%ld, m=%ld\n", pid, n, m);
		matrix[pid].resize(n);
		for (int i = 0; i < n; i++) {
			matrix[pid][i].resize(m);
			in.read((char*)matrix[pid][i].data(), sizeof(int) * m);
		}
	}
	in.close();
	return matrix;
}

void align_all_repeats(const ZigOptions &opt,
	const string &name_t, vector<TandemRepeats> &lr_t, int t_len, const char *t_seq,
	const string &name_q, vector<TandemRepeats> &lr_q, int q_len, const char *q_seq)
{
	int n = 0, m = 0;
	int t_sum = 0, q_sum = 0;
	for (TandemRepeats &lr: lr_t) {
		lr.base_idx = n;
		n += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			t_sum += r.end - r.beg;
		}
	}
	for (TandemRepeats &lr: lr_q) {
		lr.base_idx = m;
		m += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			q_sum += r.end - r.beg;
		}
	}
	fprintf(stderr, "T: %d repeats, fraction=%.2f %%\n", n, 100.0 * t_sum / t_len);
	fprintf(stderr, "Q: %d repeats, fraction=%.2f %%\n", m, 100.0 * q_sum / q_len);

	// Build score matrix of repeat units
	double t_cpu = cputime(), t_real = realtime();
	vector<vector<int>> matrix(n); // FIXME: memory inefficient for this sparse matrix
	for (int i = 0; i < n; i++) matrix[i].resize(m, -INF);
	for (const TandemRepeats &tr_t : lr_t) {
		int idx_t = tr_t.base_idx;
		int pid_t = tr_t.pid;
		for (const TandemRepeats &tr_q : lr_q) {
			int idx_q = tr_q.base_idx;
			int pid_q = tr_q.pid;
			if (pid_q != pid_t) continue;

			#pragma omp parallel for
			for (int i = 0; i < tr_t.repeats.size(); i++) {
				const RepInterval &t = tr_t.repeats[i];
				for (int j = 0; j < tr_q.repeats.size(); j++) {
					const RepInterval &q = tr_q.repeats[j];
					matrix[idx_t + i][idx_q + j] = global_alignment(
						opt.unit_aln, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
				}
			}
		}
	}
	fprintf(stderr, "Build scoring matrix: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);

	// TODO: set the penalty of deleting units
	// sqrt(2 * L * max_div * mis_penalty)
	t_cpu = cputime(); t_real = realtime();
	const int DEL_UNIT = opt.del_unit_pen;
	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;
	vector<vector<int>> dp(n + 1);
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) {
		dp[i].resize(m + 1, -INF);
		bt[i].resize(m + 1, 0);
	}
	dp[0][0] = 0;
	for (int j = 1; j <= m; j++) {
		dp[0][j] = DEL_UNIT * j;
	}
	for (int i = 1; i <= n; i++) {
		dp[i][0] = DEL_UNIT * i;
		for (int j = 1; j <= m; j++) {
			int v = dp[i-1][j] + DEL_UNIT;
			int h = dp[i][j-1] + DEL_UNIT;
			int d = dp[i-1][j-1] + matrix[i-1][j-1];
			if (v > dp[i][j]) {
				dp[i][j] = v;
				bt[i][j] = VERTICAL;
			}
			if (h > dp[i][j]) {
				dp[i][j] = h;
				bt[i][j] = HORIZONTAL;
			}
			if (d > dp[i][j]) {
				dp[i][j] = d;
				bt[i][j] = DIAGONAL;
			}
		}
	}
	fprintf(stderr, "DP: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);

	// Construct alignment topology
	t_cpu = cputime(); t_real = realtime();
	vector<pair<int,int>> aln;
	vector<int> t_del;
	vector<int> q_del;
	int pi = n, pj = m;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case DIAGONAL:
			pi--;
			pj--;
			aln.emplace_back(make_pair(pi, pj));
			break;
		case VERTICAL:
			pi--;
			t_del.push_back(pi);
			break;
		case HORIZONTAL:
			pj--;
			q_del.push_back(pj);
			break;
		default:
			break;
		}
	}
	reverse(aln.begin(), aln.end());
	reverse(t_del.begin(), t_del.end());
	reverse(q_del.begin(), q_del.end());

	int aln_len = 0, t_del_len = 0, q_del_len = 0;
	FILE *fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_aln.tsv").c_str(), "w") :nullptr;
	if (fo) {
		fprintf(fo, "%s\t%s\t%s\t%s\t", "t_id", "t_beg", "t_end", "t_len");
		fprintf(fo, "%s\t%s\t%s\t%s\t", "q_id", "q_beg", "q_end", "q_len");
		fprintf(fo, "%s\t%s\t%s\n", "match", "mismatch", "gap");
	}
	vector<RepInterval> agg_t;
	for (const TandemRepeats &lr: lr_t) {
		agg_t.insert(agg_t.end(), lr.repeats.begin(), lr.repeats.end());
	}
	vector<RepInterval> agg_q;
	for (const TandemRepeats &lr: lr_q) {
		agg_q.insert(agg_q.end(), lr.repeats.begin(), lr.repeats.end());
	}

	vector<int> rep_id_map(t_len + 1, -1);
	for (const auto &lr :lr_t) {
		for (const auto &r: lr.repeats) {
			rep_id_map[r.beg] = lr.pid; // Mark the position one base before the repeat
		}
	}

	AlnSta mut; // Small mutations
	for (const pair<int,int> &p: aln) {
		int i = p.first, j = p.second;
		const RepInterval &t = agg_t[i];
		const RepInterval &q = agg_q[j];
		aln_len += min(t.end - t.beg, q.end - q.beg);
		AlnSta res = global_cigar(opt, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
		mut += res;
		if (fo) {
			fprintf(fo, "%d\t%d\t%d\t%d\t", i, t.beg, t.end, t.end - t.beg);
			fprintf(fo, "%d\t%d\t%d\t%d\t", j, q.beg, q.end, q.end - q.beg);
			fprintf(fo, "%d\t%d\t%d\n", res.match, res.mismatch, res.ins + res.del);
		}
	}
	if (fo) fclose(fo);

	double mut_rate = 100.0 * (mut.mismatch + mut.ins + mut.del) / aln_len;
	fprintf(stderr, "Alignment between repeats:\n");
	fprintf(stderr, "    aligned length: %d, map ratio: %.2f %%\n", aln_len, 100.0 * aln_len / min(t_len, q_len));
	fprintf(stderr, "    matches: %d, mismatches: %d, insertions: %d, deletions: %d, mutation rate: %.2f %%\n",
		mut.match, mut.mismatch, mut.ins, mut.del, mut_rate);

	// Target deletion
	fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_t_del.tsv").c_str(), "w") :nullptr;
	if (fo) fprintf(fo, "%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len");
	for (int i: t_del) {
		const RepInterval &t = agg_t[i];
		if (fo) fprintf(fo, "%d\t%d\t%d\t%d\n", i, t.beg, t.end, t.end - t.beg);
		t_del_len += t.end - t.beg;
	}
	if (fo) fclose(fo);
	fprintf(stderr, "    target deletion length: %d, percentage: %.2f %%\n", t_del_len, 100.0 * t_del_len / t_len);

	// Query deletion
	fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_q_del.tsv").c_str(), "w") :nullptr;
	if (fo) fprintf(fo, "%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len");
	for (int j: q_del) {
		const RepInterval &q = agg_q[j];
		if (fo) fprintf(fo, "%d\t%d\t%d\t%d\n", j, q.beg, q.end, q.end - q.beg);
		q_del_len += q.end - q.beg;
	}
	if (fo) fclose(fo);
	fprintf(stderr, "    query deletion length: %d, percentage: %.2f %%\n", q_del_len, 100.0 * q_del_len / q_len);

	fprintf(stderr, "Construct topology: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);

	// Finalizing the alignment
	vector<RepInterval> aln_t, aln_q;
	vector<RepInterval> del_t, del_q;
	for (auto &pair: aln) {
		aln_t.push_back(agg_t[pair.first]);
		aln_q.push_back(agg_q[pair.second]);
	}
	for (int i: t_del) {
		del_t.push_back(agg_t[i]);
	}
	for (int i: q_del) {
		del_q.push_back(agg_q[i]);
	}

	{
		for (int i = 1; i < aln_t.size(); i++) {
			assert(aln_t[i].beg >= aln_t[i-1].end);
			assert(aln_q[i].beg >= aln_q[i-1].end);
		}
		for (int i = 1; i < del_t.size(); i++) {
			assert(del_t[i].end >= del_t[i-1].beg);
		}
		for (int i = 1; i < del_q.size(); i++) {
			assert(del_q[i].end >= del_q[i-1].beg);
		}
	}

	string final_ct, final_cq;
	int jt = 0, jq = 0;
	int long_gap_cnt = 0;
	for (int i = 0; i <= aln_q.size(); i++) {
		int t_beg = i == 0 ?0 :aln_t[i-1].end, t_end = i == aln_q.size() ?t_len :aln_t[i].beg;
		int q_beg = i == 0 ?0 :aln_q[i-1].end, q_end = i == aln_q.size() ?q_len :aln_q[i].beg;
		int gap_t_len = t_end - t_beg, gap_q_len = q_end - q_beg;
		// fprintf(stderr, "[%d,%d) t_len=%d -> [%d,%d) q_len=%d\n", t_beg, t_end, len_t, q_beg, q_end, len_q);

		// Find deleted segments within the gap
		vector<RepInterval> sub_t;
		for (; jt < del_t.size(); jt++) {
			int b = del_t[jt].beg, e = del_t[jt].end;
			if (b >= t_beg and e <= t_end) {
				RepInterval x = del_t[jt];
				x.beg -= t_beg;
				x.end -= t_beg;
				sub_t.push_back(x);
			} else if (b >= t_end) {
				break;
			}
		}
		vector<RepInterval> sub_q;
		for (; jq < del_q.size(); jq++) {
			int b = del_q[jq].beg, e = del_q[jq].end;
			if (b >= q_beg and e <= q_end) {
				RepInterval x = del_q[jq];
				x.beg -= q_beg;
				x.end -= q_beg;
				sub_q.push_back(x);
			} else if (b >= q_end) {
				break;
			}
		}
		// Caution: deletion-insertion runs might be caused by inaccurate partitioning

		// Collect the left segments
		string left_t;
		vector<int> map_t;
		const char *st = t_seq + t_beg;
		for (int j = 0; j <= sub_t.size(); j++) {
			int b = (j == 0) ?0 :sub_t[j-1].end;
			int e = (j == sub_t.size()) ?gap_t_len :sub_t[j].beg;
			for (int k = b; k < e; k++) {
				left_t.push_back(st[k]);
				map_t.push_back(k);
			}
		}
		left_t.push_back('$'); // Guard
		map_t.push_back(gap_t_len);

		string left_q;
		vector<int> map_q;
		const char *sq = q_seq + q_beg;
		for (int j = 0; j <= sub_q.size(); j++) {
			int b = (j == 0) ?0 :sub_q[j-1].end;
			int e = (j == sub_q.size()) ?gap_q_len :sub_q[j].beg;
			for (int k = b; k < e; k++) {
				left_q.push_back(sq[k]);
				map_q.push_back(k);
			}
		}
		left_q.push_back('$');
		map_q.push_back(gap_q_len);

		// Be careful with long gaps
		if (left_t.length() > 5000 or left_q.length() > 5000) {
			long_gap_cnt++;
		}

		// FIXME: use SI scoring matrix
		AlnSta as = global_cigar(opt, left_t.length(), left_t.data(), left_q.length(), left_q.data());
		const string &ext_t = as.ext_a;
		const string &ext_q = as.ext_b;

		int last_t = 0, last_q = 0;
		int pnt_t = 0, pnt_q = 0;
		assert(final_ct.size() == final_cq.size());
		int kt = 0, kq = 0;
		for (int j = 0; j < ext_t.length(); j++) {
			if (ext_t[j] != '-') {
				assert(ext_t[j] == left_t[pnt_t]);
				int pos_t = map_t[pnt_t];
				if (ext_t[j] != '$') assert(ext_t[j] == st[pos_t]);
				// Duplication deletions between bases
				for (; kt < sub_t.size(); kt++) {
					int b = sub_t[kt].beg, e = sub_t[kt].end;
					if (b >= last_t and e <= pos_t) {
						final_ct.push_back('[');
						final_cq.push_back('[');
						for (int c = b; c < e; c++) {
							final_ct.push_back(st[c]);
							final_cq.push_back('-');
						}
						final_ct.push_back(']');
						final_cq.push_back(']');
					} else if (b >= pos_t) break;
				}
				last_t = pos_t;
				pnt_t++;
			}
			if (ext_q[j] != '-') {
				assert(ext_q[j] == left_q[pnt_q]);
				int pos_q = map_q[pnt_q];
				if (ext_q[j] != '$') assert(ext_q[j] == sq[pos_q]);
				for (; kq < sub_q.size(); kq++) {
					int b = sub_q[kq].beg, e = sub_q[kq].end;
					if (b >= last_q and e <= pos_q) {
						final_ct.push_back('[');
						final_cq.push_back('[');
						for (int c = b; c < e; c++) {
							final_ct.push_back('-');
							final_cq.push_back(sq[c]);
						}
						final_ct.push_back(']');
						final_cq.push_back(']');
					} else if (b >= pos_q) break;
				}
				last_q = pos_q;
				pnt_q++;
			}
			if (ext_t[j] != '$') final_ct.push_back(ext_t[j]);
			if (ext_q[j] != '$') final_cq.push_back(ext_q[j]);
		}
		assert(final_cq.size() == final_ct.size());

		if (i < aln_q.size()) {
			as = global_cigar(opt, aln_t[i].end - aln_t[i].beg, t_seq + aln_t[i].beg,
				aln_q[i].end - aln_q[i].beg, q_seq + aln_q[i].beg);
			final_ct.push_back('[');
			final_ct += as.ext_a;
			final_ct.push_back(']');
			final_cq.push_back('[');
			final_cq += as.ext_b;
			final_cq.push_back(']');
		}
	}
	fprintf(stderr, "Found %d long gaps\n", long_gap_cnt);

	// Sanity check
	int i = 0;
	for (char c: final_ct) {
		if (c == '[' or c == ']' or c == '-') {
			continue;
		}
		assert(c == t_seq[i++]);
	}
	assert(i == t_len);

	i = 0;
	for (char c: final_cq) {
		if (c == '[' or c == ']' or c == '-') {
			continue;
		}
		assert(c == q_seq[i++]);
	}
	assert(i == q_len);

	// fprintf(stdout, "%s\n", final_ct.data());
	// fprintf(stdout, "%s\n", final_cq.data());

	int cnt_mat = 0, cnt_mis = 0, cnt_del = 0, cnt_ins = 0;
	for (i = 0; i < final_ct.length(); i++) {
		if (final_ct[i] == '-') cnt_ins++;
		else if (final_cq[i] == '-') cnt_del++;
		else if (final_cq[i] == final_ct[i]) cnt_mat++;
		else cnt_mis++;
	}

	fprintf(stderr, "Finalized alignment: ");
	int cnt_aln = cnt_mis + cnt_mat;
	fprintf(stderr, "    aligned length: %d, map ratio: %.2f %%\n", cnt_aln, 100.0 * cnt_aln / min(t_len, q_len));
	fprintf(stderr, "    mismatch: %d, percentage: %.2f %%\n", cnt_mis, 100.0 * cnt_mis / cnt_aln);
	fprintf(stderr, "    deletion: %d, percentage: %.2f %%\n", cnt_del, 100.0 * cnt_del / t_len);
	fprintf(stderr, "    insertion: %d, percentage: %.2f %%\n", cnt_ins, 100.0 * cnt_ins / q_len);
	fprintf(stderr, "\n");

	// Output CIGAR
	extended_paf_format(name_t, t_len, final_ct, name_q, q_len, final_cq, rep_id_map);
}

void align_etr(const ZigOptions &opt, const char *fn_t, const char *fn_q)
{
	FastaSeq target = input_fasta_seq(fn_t);
	FastaSeq query = input_fasta_seq(fn_q);
	const string &name_t = target.name, &name_q = query.name;
	int len_t = target.seq.length(), len_q = query.seq.length();
	const char *seq_t = target.seq.data(), *seq_q = query.seq.data();

	vector<TandemRepeats> lr_t, lr_q;
	int n_class;
	if (opt.dev_label[0] == '1') {
		lr_t = collect_long_repeats(opt, len_t, seq_t);
		lr_q = collect_long_repeats(opt, len_q, seq_q);
		n_class = categorize_repeats(opt, lr_t, seq_t, lr_q, seq_q);
		n_class = normalize_repeats(opt, n_class, len_t, seq_t, lr_t, len_q, seq_q, lr_q);
		if (opt.log_prefix) {
			string reps_fn_t = string(opt.log_prefix) + "_t_reps.bin";
			string reps_fn_q = string(opt.log_prefix) + "_q_reps.bin";
			output_repeats(reps_fn_t, lr_t);
			output_repeats(reps_fn_q, lr_q);
		}
	} else {
		string reps_fn_t = string(opt.log_prefix) + "_t_reps.bin";
		string reps_fn_q = string(opt.log_prefix) + "_q_reps.bin";
		lr_t = input_repeats(reps_fn_t);
		lr_q = input_repeats(reps_fn_q);
		n_class = 0;
		for (const TandemRepeats &tr: lr_t) n_class = max(n_class, tr.pid + 1);
		for (const TandemRepeats &tr: lr_q) n_class = max(n_class, tr.pid + 1);
	}

	// print_tr_by_order(lr_t,lr_q);
	// print_tr_by_class(n_class, lr_t, lr_q);

	vector<vector<vector<int>>> matrix;
	if (opt.dev_label[1] == '1') {
		matrix = build_score_matrix(opt, n_class, lr_t, len_t, seq_t, lr_q, len_q, seq_q);
		string matrix_fn = string(opt.log_prefix) + "_matrix.bin";
		output_matrix(matrix_fn, matrix);
	} else {
		string matrix_fn = string(opt.log_prefix) + "_matrix.bin";
		matrix = input_matrix(matrix_fn);
	}
	exit(1);

	align_all_repeats(opt, name_t, lr_t, len_t, seq_t, name_q, lr_q, len_q, seq_q);
}
