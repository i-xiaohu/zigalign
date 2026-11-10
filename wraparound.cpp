#include "wraparound.h"
#include "ksw2.h"

int global_alignment(const ScoreParameters &sp, int n, const char *a, int m, const char *b)
{
    const int MAT_S = sp.mat;
    const int MIS_P = sp.mis;
    const int GAP_O = sp.gap_o;
    const int GAP_E = sp.gap_e;
    // -mat_len * mat_score ~= gap_len * gap_extension
    const double GAP_RATIO = -1.0 * MAT_S / GAP_E;
    int w = min(n, m) * GAP_RATIO;

    vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
    vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
    prev_H[0] = 0;
    for (int i = 1; i <= m; i++) {
        prev_H[i] = prev_E[i] = GAP_O + i * GAP_E;
    }

    int ret = -INF;
    for (int i = 1; i <= n; i++) {
        int beg = max(i - w, 1), end = min(i + w, m);
        if (beg > end) break;
        if (beg == 1) curr_H[beg-1] = GAP_O + i * GAP_E;
        else curr_H[beg-1] = -INF;
        int F = -INF;
        int max_score = -INF;
        for (int j = beg; j <= end; j++) {
            F = max(F, curr_H[j-1] + GAP_O) + GAP_E;
            curr_E[j] = max(prev_E[j], prev_H[j] + GAP_O) + GAP_E;
            int M = prev_H[j-1] + (a[i-1] == b[j-1] ?MAT_S :MIS_P);
            curr_H[j] = max(F, M);
            curr_H[j] = max(curr_H[j], curr_E[j]);
            max_score = max(max_score, curr_H[j]);
        }
        if (end < m) curr_H[end+1] = curr_E[end+1] = -INF;
        if (max_score == -INF) break;
        if (i == n and end == m) ret = curr_H[m];

        swap(prev_H, curr_H);
        swap(prev_E, curr_E);
    }
    return ret;
}

AlnSta global_cigar(const ScoreParameters &sp, int n, const char *a, int m, const char *b, int bw)
{
	if (1L * n * m / 1e9 > 64) {
		fprintf(stderr, "[%s] Memery requirement > 64G !!!\n", __func__);
		exit(1);
	}

	const int MAT_S = sp.mat;
	const int MIS_P = sp.mis;
	const int GAP_O = sp.gap_o;
	const int GAP_E = sp.gap_e;

	vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
	vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
	prev_H[0] = 0;
	for (int i = 1; i <= m; i++) prev_H[i] = prev_E[i] = GAP_O + i * GAP_E;
	// TODO: reduce memory cost or use KSW
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) bt[i].resize(m + 1, 0);
	for (int j = 1; j <= m; j++) bt[0][j] = BT_HORIZONTAL;
	for (int i = 1; i <= n; i++) bt[i][0] = BT_VERTICAL;

	int beg = 1, end = m;
	for (int i = 1; i <= n; i++) {
		beg = max(1, i - bw);
		end = min(m, i + bw);
		if (beg > end) {
			// Band moves out of the matrix from the right side
			assert(end == m);
			for (int k = i; k <= n; k++) bt[k][m] = BT_VERTICAL;
			break;
		}
		if (beg == 1) curr_H[beg-1] = GAP_O + i * GAP_E;
		else curr_H[beg-1] = -INF;
		int F = -INF;
		for (int j = beg; j <= end; j++) {
			F = max(F, curr_H[j-1] + GAP_O) + GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + GAP_O) + GAP_E;
			if (F > curr_E[j]) {
				curr_H[j] = F;
				bt[i][j] = BT_HORIZONTAL;
			} else {
				curr_H[j] = curr_E[j];
				bt[i][j] = BT_VERTICAL;
			}
			if (a[i-1] != UNIT_SI_CHAR and b[j-1] != UNIT_SI_CHAR) {
				int M = prev_H[j-1] + (a[i-1] == b[j-1] ?MAT_S :MIS_P);
				if (M > curr_H[j]) {
					curr_H[j] = M;
					bt[i][j] = BT_DIAGONAL;
				}
			}
		}
		if (end < m) curr_H[end+1] = curr_E[end+1] = -INF;

		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}
	if (end < m) {
		// Band didn't reach the bottom-right corner
		for (int j = end+1; j <= m; j++) bt[n][j] = BT_HORIZONTAL;
	}

	AlnSta ret;
	int pi = n, pj = m;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case BT_DIAGONAL:
			pi--;
			pj--;
			if (a[pi] == b[pj]) ret.mat++;
			else ret.mis++;
			ret.ext_a += a[pi];
			ret.ext_b += b[pj];
			break;
		case BT_VERTICAL:
			pi--;
			ret.del++;
			ret.ext_a += a[pi];
			ret.ext_b += '-';
			break;
		case BT_HORIZONTAL:
			pj--;
			ret.ins++;
			ret.ext_a += '-';
			ret.ext_b += b[pj];
			break;
		default:
			break;
		}
	}
	assert(pi == 0 and pj == 0);
	reverse(ret.ext_a.begin(), ret.ext_a.end());
	reverse(ret.ext_b.begin(), ret.ext_b.end());

	// Sanity check
	{
		assert(ret.ext_a.length() == ret.ext_b.length());
		int i = 0;
		for (char c: ret.ext_a) {
			if (c != '-') {
				assert(c == a[i++]);
			}
		}
		assert(i == n);

		i = 0;
		for (char c: ret.ext_b) {
			if (c != '-') {
				assert(c == b[i++]);
			}
		}
		assert(i == m);

		for (i = 0; i < ret.ext_a.length(); i++) {
			if (ret.ext_a[i] == UNIT_SI_CHAR) assert(ret.ext_b[i] == '-');
			if (ret.ext_b[i] == UNIT_SI_CHAR) assert(ret.ext_a[i] == '-');
		}
	}
	return ret;
}

AlnSta ksw_wrapper(const ScoreParameters &sp, int n, const char *a, int m, const char *b, int bw) {
	if (n == 0 or m == 0) {
		AlnSta ret;
		if (n == 0) {
			for (int i = 0; i < m; i++) {
				ret.ext_a += '-';
				ret.ext_b += b[i];
			}
		}
		if (m == 0) {
			for (int i = 0; i < n; i++) {
				ret.ext_a += a[i];
				ret.ext_b += '-';
			}
		}
		return ret;
	}

	const int MAT_S = abs(sp.mat);
	const int MIS_P = -abs(sp.mis);
	const int GAP_O = abs(sp.gap_o);
	const int GAP_E = abs(sp.gap_e);

	int8_t mat[16];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i * 4 + j] = (i == j) ?MAT_S :MIS_P;
		}
	}

	int8_t base2num[256];
	memset(base2num, 0, sizeof(base2num));
	base2num['A'] = base2num['a'] = 0;
	base2num['C'] = base2num['c'] = 1;
	base2num['G'] = base2num['g'] = 2;
	base2num['T'] = base2num['t'] = 3;

	int tl = n, ql = m;
	uint8_t *ts = (uint8_t*) malloc(tl);
	uint8_t *qs = (uint8_t*) malloc(ql);
	for (int i = 0; i < tl; i++) ts[i] = base2num[(uint8_t)a[i]];
	for (int i = 0; i < ql; i++) qs[i] = base2num[(uint8_t)b[i]];

	ksw_extz_t ez;
	memset(&ez, 0, sizeof(ksw_extz_t));
	// extz reported corrupted memory
	// ksw_extz(0, ql, qs, tl, ts, 4, mat, GAP_O, GAP_E, bw, -1, 0, &ez);
	// ksw_gg reported segmental fault
	// ksw_gg(0, ql, qs, tl, ts,4, mat,GAP_O, GAP_E, bw, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
#ifdef __SSE2__
	ksw_extz2_sse(0, ql, qs, tl, ts, 4, mat, GAP_O, GAP_E, bw, -1, 0, -1, &ez);
#else
	ksw_gg2(0, ql, qs, tl, ts,4, mat,GAP_O, GAP_E, bw, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
#endif

	AlnSta ret;
	int pi = 0, pj = 0;
	for (int i = 0; i < ez.n_cigar; i++) {
		int op = ez.cigar[i] & 0xf;
		int num = ez.cigar[i] >> 4;
		if (op == KSW_CIGAR_MATCH) {
			for (int k = 0; k < num; k++) {
				ret.ext_a += a[pi++];
				ret.ext_b += b[pj++];
			}
		} else if (op == KSW_CIGAR_INS) {
			for (int k = 0; k < num; k++) {
				ret.ext_a += '-';
				ret.ext_b += b[pj++];
			}
		} else if (op == KSW_CIGAR_DEL) {
			for (int k = 0; k < num; k++) {
				ret.ext_a += a[pi++];
				ret.ext_b += '-';
			}
		}
	}

	// Sanity check
	{
		assert(pi == n and pj == m);
		assert(ret.ext_a.length() == ret.ext_b.length());
		int i = 0;
		for (char c: ret.ext_a) {
			if (c != '-') {
				assert(c == a[i++]);
			}
		}
		assert(i == n);

		i = 0;
		for (char c: ret.ext_b) {
			if (c != '-') {
				assert(c == b[i++]);
			}
		}
		assert(i == m);
	}
	kfree(0, ez.cigar);
	free(ts);
	free(qs);

	return ret;
}

struct WdpCell {
	int E, F, H;

	WdpCell() {
		E = F = H = -INF;
	}
};

vector<RepInterval> wraparound_dp(const ScoreParameters &sp, int pat_len, const char *pat_seq, int que_len, const char *que_seq)
{
	// Classical SI score matrix
	const int MAT_S = sp.mat;
	const int MIS_P = sp.mis;
	const int GAP_O = sp.gap_o;
	const int GAP_E = sp.gap_e;

	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;
	const int WRAP_1 = 4;
	const int WRAP_2 = 5;

	vector<WdpCell> prev_dp( pat_len + 1);
	vector<WdpCell> curr_dp(pat_len + 1);
	vector<vector<uint8_t>> bt(que_len + 1);
	for (int i = 0; i <= que_len; i++) {
		bt[i].resize(pat_len + 1, 0);
	}

	prev_dp[0].H = 0;
	for (int j = 1; j <= pat_len; j++) {
		prev_dp[j].H = prev_dp[j].E = GAP_O + j * GAP_E;
		bt[0][j] = HORIZONTAL;
	}
	for (int i = 1; i <= que_len; i++) {
		// First pass
		curr_dp[0].H = GAP_O + i * GAP_E;
		bt[i][0] = VERTICAL;
		curr_dp[1].F = max(prev_dp[pat_len].H + GAP_O, prev_dp[pat_len].F) + GAP_E; // Wrap
		curr_dp[1].E = max(prev_dp[1].H + GAP_O, prev_dp[1].E) + GAP_E;
		int tmp = que_seq[i-1] == pat_seq[0] ?MAT_S :MIS_P;
		int d1 = prev_dp[0].H + tmp; // Leading deletions
		int d2 = prev_dp[pat_len].H + tmp; // Wrap
		curr_dp[1].H = -INF;
		if (curr_dp[1].F > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].F;
			bt[i][1] = WRAP_1;
		}
		if (curr_dp[1].E > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].E;
			bt[i][1] = VERTICAL;
		}
		if (d1 > curr_dp[1].H) {
			curr_dp[1].H = d1;
			bt[i][1] = DIAGONAL;
		}
		if (d2 > curr_dp[1].H) {
			curr_dp[1].H = d2;
			bt[i][1] = WRAP_1;
		}
		for (int j = 2; j <= pat_len; j++) {
			curr_dp[j].F = max(curr_dp[j-1].H + GAP_O, curr_dp[j-1].F) + GAP_E;
			curr_dp[j].E = max(prev_dp[j].H + GAP_O, prev_dp[j].E) + GAP_E;
			int d = prev_dp[j-1].H + (que_seq[i-1] == pat_seq[j-1] ?MAT_S : MIS_P);
			curr_dp[j].H = -INF;
			if (curr_dp[j].F > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].F;
				bt[i][j] = HORIZONTAL;
			}
			if (curr_dp[j].E > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].E;
				bt[i][j] = VERTICAL;
			}
			if (d > curr_dp[j].H) {
				curr_dp[j].H = d;
				bt[i][j] = DIAGONAL;
			}
		}
		// Second pass
		tmp = curr_dp[pat_len].H + GAP_O + GAP_E;
		if (tmp > curr_dp[1].F) {
			curr_dp[1].F = tmp;
		}
		if (curr_dp[1].F > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].F;
			bt[i][1] = WRAP_2;
		}
		for (int j = 2; j <= pat_len; j++) {
			tmp = max(curr_dp[j-1].H + GAP_O, curr_dp[j-1].F) + GAP_E;
			curr_dp[j].F = max(tmp, curr_dp[j].F);
			if (curr_dp[j].F > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].F;
				bt[i][j] = HORIZONTAL;
			}
		}
		swap(prev_dp, curr_dp);
	}
	vector<RepInterval> segments;
	RepInterval tmp;
	tmp.end = que_len;
	int ti = que_len, tj = pat_len;
	while (ti > 0 or tj > 0) {
		switch (bt[ti][tj]) {
			case HORIZONTAL:
				tmp.gap++;
				tj--;
				break;
			case VERTICAL:
				tmp.gap++;
				ti--;
				break;
			case DIAGONAL:
				tmp.mis += (que_seq[ti-1] != pat_seq[tj-1]);
				tmp.mat += (que_seq[ti-1] == pat_seq[tj-1]);
				ti--;
				tj--;
				break;
			case WRAP_1:
				assert(tj == 1);
				tmp.mis += (que_seq[ti-1] != pat_seq[tj-1]);
				tmp.mat += (que_seq[ti-1] == pat_seq[tj-1]);
				tmp.beg = ti - 1;
				ti--;
				tj = pat_len;
				segments.push_back(tmp);
				tmp = RepInterval();
				tmp.end = ti;
				break;
			case WRAP_2:
				assert(tj == 1);
				tmp.gap++;
				tmp.beg = ti - 1;
				tj = pat_len;
				segments.push_back(tmp);
				tmp = RepInterval();
				tmp.end = ti - 1;
				break;
			default:
				break;
		}
	}
	tmp.beg = ti;
	segments.push_back(tmp);
	reverse(segments.begin(), segments.end());
	return segments;
}

// Banding is not applicable
SgResult semi_global(const ScoreParameters &sc, int n, const char *a, int m, const char *b)
{
	const int MAT_S = sc.mat;
	const int MIS_P = sc.mis;
	const int GAP_O = sc.gap_o;
	const int GAP_E = sc.gap_e;

	vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
	vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
	for (int i = 0; i <= m; i++) {
		prev_H[i] = prev_E[i] = 0;
	}
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) {
		bt[i].resize(m + 1, 0);
	}

	for (int i = 1; i <= n; i++) {
		curr_H[0] = curr_E[0] = GAP_O + i * GAP_E;
		int F = -INF;
		for (int j = 1; j <= m; j++) {
			F = max(F, curr_H[j-1] + GAP_O) + GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + GAP_O) + GAP_E;
			int M = prev_H[j-1] + (a[i-1] == b[j-1] ?MAT_S :MIS_P);
			if (F > M) {
				curr_H[j] = F;
				bt[i][j] = BT_HORIZONTAL;
			} else {
				curr_H[j] = M;
				bt[i][j] = BT_DIAGONAL;
			}
			if (curr_E[j] > curr_H[j]) {
				curr_H[j] = curr_E[j];
				bt[i][j] = BT_VERTICAL;
			}
		}
		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}

	SgResult ret;
	for (int j = 1; j <= m; j++) {
		if (prev_H[j] > ret.score) {
			ret.score = prev_H[j];
			ret.end = j; // 0-based open
		}
	}

	int pi = n, pj = ret.end;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case BT_DIAGONAL:
			pi--;
			pj--;
			ret.mis += (a[pi] != b[pj]);
			break;
		case BT_VERTICAL:
			pi--;
			ret.gap++;
			break;
		case BT_HORIZONTAL:
			pj--;
			ret.gap++;
			break;
		default:
			break;
		}
	}
	ret.gap += pi;
	ret.beg = pj; // 0-based closed
	return ret;
}