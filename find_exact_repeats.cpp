//
// Created by ixiaohu on 2026/1/16.
//

#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <ctime>
using namespace std;

const int MIN_UNIT_LEN = 3;

// Longest prefix-substring match within S
// Let K = LPS[i], S[0:K) = S[i:i+K)
vector<int> calc_lps(const string &s) {
	const int n = s.length();
	vector<int> lps(n); // lps[i] = max{len | S[0:len)=S[i:i+len)}
	lps[0] = 0;
	int len = 0;
	while (1+len < n and s[0+len] == s[1+len]) len++;
	lps[1] = len;

	int k = 1;
	for (int i = 2; i < n; i++) {
		int u = k + lps[k];
		if (u > i) {
			int t = u - i;
			// S[i:u) equals to S[l-t:l), which serves as a bridge
			for (int j = 0; j < t; j++) {
				assert(s[i+j] == s[lps[k]-t+j]);
			}
			// Compare with prefix
			int x = lps[k] - t;
			if (lps[x] < t) {
				lps[i] = lps[x];
			} else {
				// The match can further extend
				len = t;
				while (i + len < n and s[0+len] == s[i+len]) len++;
				lps[i] = len;
			}
		} else {
			len = 0;
			while (i + len < n and s[0+len] == s[i+len]) len++;
			lps[i] = len;
		}
		if (i + lps[i] > k + lps[k]) {
			k = i;
		}
	}

	// Check correctness
	for (int i = 1; i < n; i++) {
//		printf("i=%d, lps=%d\n", i, lps[i]);
		len = lps[i];
		for (int j = 0; j < len; j++) {
			assert(s[0+j] == s[i+j]);
		}
		assert(i+len == n or s[len] != s[i+len]);
	}
	return lps;
}

// Longest substring-prefix match between T and Q
// Let LSP[i] = K, T[i:i+K) = Q[0:K)
vector<int> calc_lsp(const string &t, const string &q) {
	const int n = t.length();
	const int m = q.length();
	vector<int> lsp(n);
	int len = 0;
	while (len < n and len < m and q[len] == t[len]) len++;
	lsp[0] = len;

	vector<int> tmp = calc_lps(q);
	int k = 0;
	for (int i = 1; i < n; i++) {
		int u = k + lsp[k];
		if (u > i) {
			// Similarly, T[i:i+t) = Q[l-t:l) is the bridge
			int tt = u - i;
			for (int j = 0; j < tt; j++) {
				assert(t[i+j] == q[lsp[k]-tt+j]);
			}
			int x = lsp[k] - tt;
			if (tmp[x] < tt) {
				lsp[i] = tmp[x];
			} else {
				len = tt;
				while (i + len < n and len < m and q[0+len] == t[i+len]) len++;
				lsp[i] = len;
			}
		} else {
			len = 0;
			while (i + len < n and len < m and q[0+len] == t[i+len]) len++;
			lsp[i] = len;
		}
		if (i + lsp[i] > k + lsp[k]) {
			k = i;
		}
	}
	// Check correctness
	for (int i = 0; i < n; i++) {
//		printf("i=%d, lsp=%d\n", i, lsp[i]);
		len = lsp[i];
		for (int j = 0; j < len; j++) {
			assert(q[0+j] == t[i+j]);
		}
		assert(i+len == n or len == m or q[len] != t[i+len]);
	}
	return lsp;
}

struct NewRep {
	int s, e, len;
	NewRep(): s(0), e(0), len(0) {}
	bool operator < (const NewRep &r) const {
		return len != r.len ?len > r.len :s < r.s;
	}
};

// Search the repetitions formed by concatenating Q and T and centered at somewhere in T
vector<NewRep> find_center_at_right(string &q, string &t)
{
	const int n = t.length();
	const int m = q.length();
	// The prefix-substring match within T
	vector<int> lps = calc_lps(t);
	// The suffix-substring match between T and Q
	reverse(t.begin(), t.end());
	reverse(q.begin(), q.end());
	vector<int> lss = calc_lsp(t, q);
	reverse(t.begin(), t.end());
	reverse(q.begin(), q.end());
	reverse(lss.begin(), lss.end());
	// Check correctness
	for (int i = 0; i < n; i++) {
		int len = lss[i];
		for (int j = 0; j < len; j++) {
			assert(q[m-1-j] == t[i-j]);
		}
	}
	vector<NewRep> ret;
	int max_len = min((n + m) / 2, n); // Maximum repeat unit length
	for (int u = MIN_UNIT_LEN; u <= max_len; u++) { // Iterate over unit length
		// The endpoints of repeats with a same length are consecutive
		// If repeat ends at i in T, then
		//     T[0:i-u] = T[u:i], i.e., LPS[u] >= i-u+1
		//     Q[m-2u+i+1:m-1] = T[i-u+1:u-1], i.e., LSS[u-1] >= 2u-i-1
		int end = min(2 * u - 2, u + (u < n ?lps[u] :0) - 1);
		int beg = max(u - 1, 2 * u - lss[u-1] - 1);
		if (beg <= end) {
			NewRep rep;
			for (int i = beg; i <= end; i++) {
				rep.s = m-2*u+i+1;
				rep.e = i;
				rep.len = u;
				ret.push_back(rep);
			}
//			fprintf(stderr, "Find %d repeats length of %d, ending at [%d,%d]\n", end-beg+1, u, beg, end);
			for (int i = beg; i <= end; i++) {
				for (int j = 0; j < i-u+1; j++) {
					assert(t[j] == t[u+j]);
				}
				for (int j = 0; j < 2*u-i-1; j++) {
					assert(q[m-2*u+i+1 + j] == t[i-u+1 + j]);
				}
			}
		}
	}
	return ret;
}

void divide_and_conquer(const string &s, int l, int r, vector<NewRep> &ans)
{
//	fprintf(stderr, "[%d,%d)\n", l, r);
	if (r - l < MIN_UNIT_LEN * 2) {
		return ;
	}
	int m = (l + r) / 2;
	string a = s.substr(l, m-l);
	string b = s.substr(m, r-m);
	vector<NewRep> r_rep = find_center_at_right(a, b);
	for (NewRep &x : r_rep) {
		x.s += l;
		x.e += m;
		assert(x.e - x.s + 1 == 2 * x.len);
		for (int i = 0; i < x.len; i++) {
			assert(s[x.s + i] == s[x.s + x.len + i]);
		}
		ans.push_back(x);
	}

	reverse(a.begin(), a.end());
	reverse(b.begin(), b.end());
	vector<NewRep> l_rep = find_center_at_right(b, a);
	for (NewRep &x : l_rep) {
		int tmp_s = a.length() - x.e - 1 + l;
		int tmp_e = b.length() - x.s - 1 + m;
		x.s = tmp_s;
		x.e = tmp_e;
		assert(x.e - x.s + 1 == 2 * x.len);
		for (int i = 0; i < x.len; i++) {
			assert(s[x.s + i] == s[x.s + x.len + i]);
		}
		// To avoid the redundant result when two strings are same
		if (x.s != l or x.len != a.length()) {
			ans.push_back(x);
		}
	}

	divide_and_conquer(s, m, r, ans);
	divide_and_conquer(s, l, m, ans);
}

int main(int argc, char *argv[]) {
	srand48(time(nullptr));
	int n = 100000;
	string s;
	for (int i = 0; i < n; i++) {
		s += "ACGT"[rand() % 4];
	}
	vector<NewRep> reps;
	divide_and_conquer(s, 0, s.length(), reps);
	sort(reps.begin(), reps.end());
	fprintf(stderr, "Found %ld perfect repeats\n", reps.size());
	for (int i = 0; i < reps.size(); i++) {
		const NewRep &r = reps[i];
		fprintf(stderr, "i=%d, length=%d, start=%d\n", i+1, r.len, r.s);
		string a = s.substr(r.s, r.len);
		string b = s.substr(r.s + r.len, r.len);
		assert(a == b);
	}
	return 0;
}