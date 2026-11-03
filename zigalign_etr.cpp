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
bool merge_tr_test(const ZigOptions &opt, const TandemRepeats &tr_a, const char *seq_a, const TandemRepeats &tr_b, const char *seq_b)
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
			SgResult sg = semi_global(opt.si, len_a, seq_a + a.beg, sb.length(), sb.data());
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
					bool to_merge = merge_tr_test(opt, t, seq, g, seq);
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
		bool to_merge = merge_tr_test(opt, t, seq_t, q, seq_q) and
						merge_tr_test(opt, q, seq_q, t, seq_t);
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

typedef vector<vector<vector<int>>> ScoreMatrix;

ScoreMatrix build_score_matrix(const ZigOptions &opt, int n_class,
	vector<TandemRepeats> &lr_t, int t_len, const char *t_seq,
	vector<TandemRepeats> &lr_q, int q_len, const char *q_seq)
{
	int n = 0, m = 0;
	int t_sum = 0, q_sum = 0;
	vector<vector<RepInterval>> reps_t(n_class);
	vector<vector<RepInterval>> reps_q(n_class);
	for (TandemRepeats &lr: lr_t) {
		n += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			reps_t[lr.pid].push_back(r);
			t_sum += r.end - r.beg;
		}
	}
	for (TandemRepeats &lr: lr_q) {
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
	ScoreMatrix matrix(n_class);
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

void output_matrix(const string &fn, const ScoreMatrix &matrix)
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
		for (int i = 0; i < n; i++) {
			assert(matrix[pid][i].size() == m);
			out.write((char*)matrix[pid][i].data(), sizeof(int) * m);
		}
	}
	out.close();
}

ScoreMatrix input_matrix(const string &fn)
{
	ifstream in(fn, ios::binary);
	assert(in.is_open());
	ScoreMatrix matrix;
	size_t n_class;
	in.read((char*)&n_class, sizeof(int64_t));
	matrix.resize(n_class);
	for (int pid = 0; pid < n_class; pid++) {
		size_t n, m;
		in.read((char*)&n, sizeof(int64_t));
		in.read((char*)&m, sizeof(int64_t));
		if (n == 0 or m == 0) continue;
		matrix[pid].resize(n);
		for (int i = 0; i < n; i++) {
			matrix[pid][i].resize(m);
			in.read((char*)matrix[pid][i].data(), sizeof(int) * m);
		}
	}
	in.close();
	return matrix;
}

struct Dp2Profile {
	int pid, idx;
	int beg, end;
	int aln; // -1 for deletion
	int del_pen; // Penalty of deletion
	Dp2Profile(int p, int i, int b, int e): pid(p), idx(i), beg(b), end(e) {
		aln = -1;
		del_pen = 0;
	}
	Dp2Profile() = default;
};

vector<Dp2Profile> construct_profile(const ZigOptions &opt, int n_class, const vector<TandemRepeats> &lr)
{
	vector<Dp2Profile> ret;
	vector<int> cnt(n_class, 0);
	for (const TandemRepeats &tr: lr) {
		for (const RepInterval &r: tr.repeats) {
			// Synchronize with the score matrix
			Dp2Profile p(tr.pid, cnt[tr.pid]++, r.beg, r.end);
			// Deletion penalty = SI_gap_o + sqrt(unit_len) * SI_gap_e
			p.del_pen = opt.si.gap_o + opt.si.gap_e * round(sqrt(r.end - r.beg));
			ret.push_back(p);
		}
	}
	for (int i = 1; i < ret.size(); i++) {
		assert(ret[i].end >= ret[i-1].beg);
	}
	return ret;
}

// This function constructs an alignment backbone for tandem repeats.
// Compared with SI+DSI profile, it prioritizes DSI region more than SI region.
void alignment_backbone(const ZigOptions &opt, const ScoreMatrix &matrix, vector<Dp2Profile> &pro_t, vector<Dp2Profile> &pro_q)
{
	double t_cpu = cputime(), t_real = realtime();
	int n = pro_t.size(), m = pro_q.size();
	vector<int> prev_dp(m + 1), curr_dp(m + 1);
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) bt[i].resize(m + 1, 0);
	prev_dp[0] = 0;
	for (int j = 1; j <= m; j++) {
		prev_dp[j] = prev_dp[j-1] + pro_q[j-1].del_pen;
		bt[0][j] = BT_HORIZONTAL;
	}
	const int MAX_UNITS_DEL = opt.max_units_del;
	int beg = 1, end = m;
	for (int i = 1; i <= n; i++) {
		const Dp2Profile &t = pro_t[i-1];
		curr_dp[0] = prev_dp[0] + t.del_pen; // The first column
		bt[i][0] = BT_VERTICAL;
		beg = max(i - MAX_UNITS_DEL, 1);
		end = min(i + MAX_UNITS_DEL, m);
		if (beg > end) {
			for (int k = i; k <= n; k++) bt[k][m] = BT_VERTICAL;
			break;
		}
		if (beg > 1) curr_dp[beg-1] = -INF; // Boundary of the band
		for (int j = beg; j <= end; j++) {
			const Dp2Profile &q = pro_q[j-1];
			int v = prev_dp[j] + t.del_pen;
			int h = curr_dp[j-1] + q.del_pen;
			int d = -INF;
			if (t.pid == q.pid) {
				assert(t.pid < matrix.size());
				assert(t.idx < matrix[t.pid].size());
				assert(q.idx < matrix[t.pid][t.idx].size());
				d = prev_dp[j-1] + matrix[t.pid][t.idx][q.idx];
			}

			curr_dp[j] = -INF; // Initialization
			if (v > curr_dp[j]) {
				curr_dp[j] = v;
				bt[i][j] = BT_VERTICAL;
			}
			if (h > curr_dp[j]) {
				curr_dp[j] = h;
				bt[i][j] = BT_HORIZONTAL;
			}
			if (d > curr_dp[j]) {
				curr_dp[j] = d;
				bt[i][j] = BT_DIAGONAL;
			}
		}
		if (end < m) curr_dp[end+1] = -INF; // Used by the next row
		swap(curr_dp, prev_dp);
	}
	if (end < m) for (int j = end+1; j <= m; j++) bt[n][j] = BT_HORIZONTAL;
	fprintf(stderr, "Construct alignment backbone in %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);
	fprintf(stderr, "  Score matrix size is (%d X %d)\n", n, m);

	int pi = n, pj = m;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case BT_DIAGONAL:
			pi--;
			pj--;
			pro_t[pi].aln = pj;
			pro_q[pj].aln = pi;
			break;
		case BT_VERTICAL:
			pi--;
			pro_t[pi].aln = -1;
			break;
		case BT_HORIZONTAL:
			pj--;
			pro_q[pj].aln = -1;
			break;
		default:
			break;
		}
	}
	assert(pi == 0 and pj == 0);

	int mat_n = 0, mat_len = 0;
	int del_n = 0, del_len = 0;
	int ins_n = 0, ins_len = 0;
	for (int i = 0; i < pro_t.size(); i++) {
		const Dp2Profile &t = pro_t[i];
		if (t.aln == -1) {
			del_n++;
			del_len += t.end - t.beg;
		} else {
			const Dp2Profile &q = pro_q[t.aln];
			assert(q.aln == i);
			mat_n++;
			mat_len += min(t.end - t.beg, q.end - q.beg);
		}
	}
	int mat2_n = 0, mat2_len = 0;
	for (int i = 0; i < pro_q.size(); i++) {
		const Dp2Profile &q = pro_q[i];
		if (q.aln == -1) {
			ins_n++;
			ins_len += q.end - q.beg;
		} else {
			const Dp2Profile &t = pro_t[q.aln];
			assert(t.aln == i);
			mat2_n++;
			mat2_len += min(t.end - t.beg, q.end - q.beg);
		}
	}
	assert(mat_len == mat2_len);
	assert(mat_n == mat2_n);

	fprintf(stderr, "  Aligned %d units, delete %d units, insert %d units\n", mat_n, del_n, ins_n);
	fprintf(stderr, "  Aligned %d bases, delete %d bases, insert %d bases\n", mat_len, del_len, ins_len);
	fprintf(stderr, "\n");
}

void output_profile(const string &fn, const vector<Dp2Profile> &pro_t, const vector<Dp2Profile> &pro_q)
{
	ofstream out(fn);
	assert(out.is_open());
	int64_t n = pro_t.size();
	out.write((char*)&n, sizeof(int64_t));
	out.write((char*) pro_t.data(), n * sizeof(Dp2Profile));
	n = pro_q.size();
	out.write((char*)&n, sizeof(int64_t));
	out.write((char*) pro_q.data(), n * sizeof(Dp2Profile));
	out.close();
}

void input_profile(const string &fn, vector<Dp2Profile> &pro_t, vector<Dp2Profile> &pro_q)
{
	ifstream in(fn);
	assert(in.is_open());
	int64_t n;
	in.read((char*)&n, sizeof(int64_t));
	pro_t.resize(n);
	in.read((char*)pro_t.data(), n * sizeof(Dp2Profile));
	in.read((char*)&n, sizeof(int64_t));
	pro_q.resize(n);
	in.read((char*)pro_q.data(), n * sizeof(Dp2Profile));
	in.close();
}

struct CigarExt
{
	int pid; // -1 for SI region
	vector<char> op;
	vector<int> num;
	CigarExt(): pid(-1) {}
};

CigarExt parse_aligned_seq(int n, const char *a, const char *b, bool neq)
{
	CigarExt ret;
	char prev_op = '@';
	int cnt = 0;
	for (int i = 0; i < n; i++) {
		char curr_op;
		if (a[i] == '-') {
			curr_op = 'I';
		} else if (b[i] == '-') {
			curr_op = 'D';
		} else {
			if (neq) {
				if (a[i] == b[i]) curr_op = '=';
				else curr_op = 'X';
			} else {
				curr_op = 'M';
			}
		}
		if (prev_op != '@' and curr_op != prev_op) {
			ret.op.push_back(prev_op);
			ret.num.push_back(cnt);
			cnt = 0;
		}
		cnt++;
		prev_op = curr_op;
	}
	if (prev_op != '@') {
		ret.op.push_back(prev_op);
		ret.num.push_back(cnt);
	}
	return ret;
}

void move_pointer(int &pi, int &pj, const CigarExt &ce) {
	for (int i = 0; i < ce.op.size(); i++) {
		if (ce.op[i] == 'I') {
			pj += ce.num[i];
		} else if (ce.op[i] == 'D') {
			pi += ce.num[i];
		} else {
			pi += ce.num[i];
			pj += ce.num[i];
		}
	}
}

// Merge SI and DSI cigar for compatibility with PAF format
CigarExt combine_cigar(const vector<CigarExt> &a)
{
	CigarExt ret;
	char prev_op = '@';
	int cnt = 0;
	for (const CigarExt &ce: a) {
		for (int i = 0; i < ce.op.size(); i++) {
			char curr_op = ce.op[i];
			if (prev_op != '@' and curr_op != prev_op) {
				ret.op.push_back(prev_op);
				ret.num.push_back(cnt);
				cnt = 0;
			}
			cnt += ce.num[i];
			prev_op = curr_op;
		}
	}
	if (prev_op != '@') {
		ret.op.push_back(prev_op);
		ret.num.push_back(cnt);
	}
	return ret;
}

string cigar_to_str(const CigarExt &ce) {
	string ret;
	for (int i = 0; i < ce.op.size(); i++) {
		ret += to_string(ce.num[i]);
		ret += ce.op[i];
	}
	return ret;
}

struct GapAndUnit {
	bool si_gap;
	string seq_t, seq_q;
	vector<Dp2Profile> del_t, del_q;
	Dp2Profile aln_t, aln_q;
	string ext_t, ext_q; // Alignment result
	GapAndUnit(bool gap): si_gap(gap), aln_t(), aln_q() {}
};

vector<CigarExt> fill_SI_gap(const ZigOptions &opt,
	int len_t, const char *seq_t, vector<Dp2Profile> &pro_t,
	int len_q, const char *seq_q, vector<Dp2Profile> &pro_q)
{
	vector<Dp2Profile> aln_t, aln_q, del_t, del_q;
	for (const Dp2Profile &p: pro_t) {
		if (p.aln == -1) del_t.push_back(p);
		else {
			aln_t.push_back(p);
			aln_q.push_back(pro_q[p.aln]);
		}
	}
	for (const Dp2Profile &p: pro_q) {
		if (p.aln == -1) del_q.push_back(p);
	}

	vector<GapAndUnit> tasks;
	int jt = 0, jq = 0;
	for (int i = 0; i <= aln_t.size(); i++) {
		int gap_t_beg = i == 0 ?0 :aln_t[i-1].end, gap_t_end = i == aln_q.size() ?len_t :aln_t[i].beg;
		int gap_q_beg = i == 0 ?0 :aln_q[i-1].end, gap_q_end = i == aln_q.size() ?len_q :aln_q[i].beg;

		// Deleted units within the gap
		GapAndUnit gap(true);
		int last_pos = gap_t_beg;
		for (; jt < del_t.size(); jt++) {
			int b = del_t[jt].beg, e = del_t[jt].end;
			if (b >= gap_t_beg and e <= gap_t_end) {
				for (int j = last_pos; j < b; j++) gap.seq_t.push_back(seq_t[j]);
				gap.seq_t.push_back(UNIT_SI_CHAR);
				Dp2Profile x = del_t[jt];
				gap.del_t.push_back(x);
				last_pos = x.end;
			} else if (b >= gap_t_end) {
				break;
			}
		}
		for (int j = last_pos; j < gap_t_end; j++) {
			gap.seq_t.push_back(seq_t[j]);
		}

		last_pos = gap_q_beg;
		for (; jq < del_q.size(); jq++) {
			int b = del_q[jq].beg, e = del_q[jq].end;
			if (b >= gap_q_beg and e <= gap_q_end) {
				for (int j = last_pos; j < b; j++) gap.seq_q.push_back(seq_q[j]);
				gap.seq_q.push_back(UNIT_SI_CHAR);
				Dp2Profile x = del_q[jq];
				gap.del_q.push_back(x);
				last_pos = x.end;
			} else if (b >= gap_q_end) {
				break;
			}
		}
		for (int j = last_pos; j < gap_q_end; j++) {
			gap.seq_q.push_back(seq_q[j]);
		}
		if (not gap.seq_t.empty() or not gap.seq_q.empty()) {
			tasks.push_back(gap);
		}

		if (i < aln_t.size()) {
			GapAndUnit unit(false);
			unit.aln_t = aln_t[i];
			unit.aln_q = aln_q[i];
			tasks.push_back(unit);
		}
	}
	assert(jt == del_t.size() and jq == del_q.size());

	const int MAX_BW = 20000;
	#pragma omp parallel for
	for (GapAndUnit &g: tasks) {
		int n, m;
		const char *t, *q;
		if (g.si_gap) {
			n = g.seq_t.length();
			m = g.seq_q.length();
			t = g.seq_t.data();
			q = g.seq_q.data();
		} else {
			n = g.aln_t.end - g.aln_t.beg;
			m = g.aln_q.end - g.aln_q.beg;
			t = seq_t + g.aln_t.beg;
			q = seq_q + g.aln_q.beg;
		}
		AlnSta as = global_cigar(opt.si, n, t, m, q, MAX_BW);
		g.ext_t = as.ext_a;
		g.ext_q = as.ext_b;
	}

	int pi = 0, pj = 0;
	vector<CigarExt> ret;
	for (const GapAndUnit &g: tasks) {
		if (g.si_gap) {
			// Insert deleted repeats back
			int last_pos = 0;
			jt = jq = 0;
			for (int i = 0; i < g.ext_t.length(); i++) {
				if (g.ext_t[i] == UNIT_SI_CHAR or g.ext_q[i] == UNIT_SI_CHAR) {
					if (i > last_pos) {
						CigarExt ce = parse_aligned_seq(i - last_pos,
							g.ext_t.data() + last_pos, g.ext_q.data() + last_pos, opt.cigar_neq);
						ce.pid = -1;
						ret.push_back(ce);
						move_pointer(pi, pj, ce);
					}
					last_pos = i + 1;
				}

				if (g.ext_t[i] == UNIT_SI_CHAR) {
					const Dp2Profile &pt = g.del_t[jt++];
					CigarExt ce;
					ce.pid = pt.pid;
					ce.op.push_back('D');
					ce.num.push_back(pt.end - pt.beg);
					ret.push_back(ce);
					assert(pi == pt.beg);
					move_pointer(pi, pj, ce);
					assert(pi == pt.end);
				} else if (g.ext_q[i] == UNIT_SI_CHAR) {
					const Dp2Profile &pq = g.del_q[jq++];
					CigarExt ce;
					ce.pid = pq.pid;
					ce.op.push_back('I');
					ce.num.push_back(pq.end - pq.beg);
					ret.push_back(ce);
					assert(pj == pq.beg);
					move_pointer(pi, pj, ce);
					assert(pj == pq.end);
				}
			}
			assert(jt == g.del_t.size() and jq == g.del_q.size());
			if (g.ext_t.length() > last_pos) {
				CigarExt ce = parse_aligned_seq(g.ext_t.length() - last_pos,
					g.ext_t.data() + last_pos, g.ext_q.data() + last_pos, opt.cigar_neq);
				ce.pid = -1;
				ret.push_back(ce);
				move_pointer(pi, pj, ce);
			}
		} else {
			CigarExt ce = parse_aligned_seq(g.ext_t.length(), g.ext_t.data(), g.ext_q.data(), opt.cigar_neq);
			ce.pid = g.aln_t.pid;
			ret.push_back(ce);
			assert(pi == g.aln_t.beg and pj == g.aln_q.beg);
			move_pointer(pi, pj, ce);
			assert(pi == g.aln_t.end and pj == g.aln_q.end);
		}
	}
	assert(pi == len_t and pj == len_q);

	// Sanity check
	{
		string t, q;
		pi = pj = 0;
		for (const CigarExt &ce: ret) {
			for (int i = 0; i < ce.op.size(); i++) {
				if (ce.op[i] == 'I') {
					for (int j = 0; j < ce.num[i]; j++) {
						q += seq_q[pj++];
					}
				} else if (ce.op[i] == 'D') {
					for (int j = 0; j < ce.num[i]; j++) {
						t += seq_t[pi++];
					}
				} else {
					for (int j = 0; j < ce.num[i]; j++) {
						t += seq_t[pi++];
						q += seq_q[pj++];
					}
				}
			}
		}
		assert(pi == len_t and pj == len_q);
		assert(string(seq_t) == t and string(seq_q) == q);
	}

	int si_mat = 0, si_mis = 0, si_ins = 0, si_del = 0;
	int u_aln = 0, u_ins = 0, u_del = 0;
	int b_mat = 0, b_mis = 0, b_ins = 0, b_del = 0;
	pi = pj = 0;
	for (const CigarExt &ce: ret) {
		if (ce.pid == -1) { // SI
			for (int i = 0; i < ce.op.size(); i++) {
				if (ce.op[i] == 'I') {
					si_ins += ce.num[i];
					pj += ce.num[i];
				} else if (ce.op[i] == 'D') {
					si_del += ce.num[i];
					pi += ce.num[i];
				} else {
					for (int j = 0; j < ce.num[i]; j++) {
						si_mat += seq_t[pi] == seq_q[pj];
						si_mis += seq_t[pi] != seq_q[pj];
						pi++;
						pj++;
					}
				}
			}
		} else { // DSI
			if (ce.op.size() == 1 and ce.op[0] == 'I') {
				u_ins++;
				b_ins += ce.num[0];
				pj += ce.num[0];
			} else if (ce.op.size() == 1 and ce.op[0] == 'D') {
				u_del++;
				b_del += ce.num[0];
				pi += ce.num[0];
			} else {
				u_aln++;
				for (int i = 0; i < ce.op.size(); i++) {
					if (ce.op[i] == 'I') {
						b_ins += ce.num[i];
						pj += ce.num[i];
					} else if (ce.op[i] == 'D') {
						b_del += ce.num[i];
						pi += ce.num[i];
					} else {
						for (int j = 0; j < ce.num[i]; j++) {
							b_mat += seq_t[pi] == seq_q[pj];
							b_mis += seq_t[pi] != seq_q[pj];
							pi++;
							pj++;
						}
					}
				}
			}
		}
	}
	assert(pi == len_t and pj == len_q);
	fprintf(stderr, "Final alignment:\n");
	fprintf(stderr, "  SI: %d mat, %d mis, %d ins, %d del\n", si_mat, si_mis, si_ins, si_del);
	fprintf(stderr, "  DSI(base): %d mat, %d mis, %d ins, %d del\n", b_mat, b_mis, b_ins, b_del);
	fprintf(stderr, "  DSI(unit): %d aln, %d ins, %d del\n", u_aln, u_ins, u_del);
	fprintf(stderr, "\n");
	return ret;
}

void print_paf_result(const char *fn, const string &t_name, int t_len, const string &q_name, int q_len, const vector<CigarExt> &cigar_ext)
{
	FILE *fo = stdout;
	if (fn) {
		fo = fopen(fn, "w");
		assert(fo != nullptr);
	}
	char strand = '+';
	int len_no_gap = 0, len_with_gap = 0;
	int mapq = 60;
	for (const CigarExt &ce: cigar_ext) {
		for (int i = 0; i < ce.op.size(); i++) {
			len_with_gap += ce.num[i];
			if (ce.op[i] != 'I' and ce.op[i] != 'D') {
				len_no_gap += ce.num[i];
			}
		}
	}
	CigarExt ce = combine_cigar(cigar_ext);
	string cigar = cigar_to_str(ce);

	fprintf(fo, "%s\t%d\t%d\t%d\t%c\t", q_name.c_str(), q_len, 0, q_len, strand);
	fprintf(fo, "%s\t%d\t%d\t%d\t", t_name.c_str(), t_len, 0, t_len);
	fprintf(fo, "%d\t%d\t%d\t", len_no_gap, len_with_gap, mapq);
	fprintf(fo, "cg:Z:%s\n", cigar.c_str());
	if (fn and fo) fclose(fo);
}

void print_dsi_result(const char *fn, const vector<CigarExt> &cigar_ext) {
	if (not fn) return;
	FILE *fo = fopen(fn, "w");
	assert(fo != nullptr);
	fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "rep_id", "t_beg", "t_end", "q_beg", "q_end", "cigar");
	int pos_t = 0, pos_q = 0;
	for (const CigarExt &ce: cigar_ext) {
		int start_t = pos_t, start_q = pos_q;
		move_pointer(pos_t, pos_q, ce);
		string cigar = cigar_to_str(ce);
		fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%s\n", ce.pid, start_t, pos_t, start_q, pos_q, cigar.c_str());
	}
	fclose(fo);
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

	ScoreMatrix matrix;
	if (opt.dev_label[1] == '1') {
		matrix = build_score_matrix(opt, n_class, lr_t, len_t, seq_t, lr_q, len_q, seq_q);
		if (opt.log_prefix) {
			string matrix_fn = string(opt.log_prefix) + "_matrix.bin";
			output_matrix(matrix_fn, matrix);
		}
	} else if (opt.dev_label[1] == '0') {
		string matrix_fn = string(opt.log_prefix) + "_matrix.bin";
		matrix = input_matrix(matrix_fn);
	}

	vector<Dp2Profile> profile_t, profile_q;
	if (opt.dev_label[2] == '1') {
		profile_t = construct_profile(opt, n_class, lr_t);
		profile_q = construct_profile(opt, n_class, lr_q);
		alignment_backbone(opt, matrix, profile_t, profile_q);
		if (opt.log_prefix) {
			string profile_fn = string(opt.log_prefix) + "_profile.bin";
			output_profile(profile_fn, profile_t, profile_q);
		}
	} else {
		string profile_fn = string(opt.log_prefix) + "_profile.bin";
		input_profile(profile_fn, profile_t, profile_q);
	}

	vector<CigarExt> cigar_ext = fill_SI_gap(opt, len_t, seq_t, profile_t, len_q, seq_q, profile_q);
	print_paf_result(opt.paf_fn, name_t, len_t, name_q, len_q, cigar_ext);
	print_dsi_result(opt.dsi_fn, cigar_ext);
}
