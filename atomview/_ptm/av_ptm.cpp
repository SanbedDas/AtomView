// C entry point for AtomView: PTM over a precomputed neighbor table, split over threads.
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <thread>
#include <vector>
#include "ptm_functions.h"

#if defined(_WIN32)
#define AV_EXPORT extern "C" __declspec(dllexport)
#else
#define AV_EXPORT extern "C" __attribute__((visibility("default")))
#endif

struct table { int k; const int64_t* idx; const double* vec; const int32_t* num; };

// idx[i*k + j]: j-th nearest neighbor of i (sorted, -1 = none); vec: minimum-image vector from i to it.
static int get_neighbours(void* vdata, size_t, size_t atom, int num, int* ordering, size_t* nbr_indices, int32_t* numbers, double (*pos)[3])
{
	const table* t = (const table*)vdata;
	pos[0][0] = pos[0][1] = pos[0][2] = 0;
	nbr_indices[0] = atom;
	numbers[0] = t->num ? t->num[atom] : 0;
	ordering[0] = 0;
	int n = 1;
	for (int j = 0; j < t->k && n < num; j++) {
		int64_t q = t->idx[atom * t->k + j];
		if (q < 0) break;
		memcpy(pos[n], &t->vec[(atom * t->k + j) * 3], 3 * sizeof(double));
		nbr_indices[n] = (size_t)q;
		numbers[n] = t->num ? t->num[q] : 0;
		ordering[n] = n;
		n++;
	}
	return n;
}

AV_EXPORT int av_ptm(int64_t n, int k, const int64_t* idx, const double* vec, const int32_t* numbers, int32_t flags,
                     int conventional, int nthreads, int32_t* type, int32_t* alloy, double* rmsd, double* scale, double* q,
                     double* F, double* inter, double* lattice, int32_t* tmpl)
{
	ptm_initialize_global();
	table t = {k, idx, vec, numbers};
	nthreads = std::max(1, std::min<int>(nthreads, (int)std::max<int64_t>(1, n / 1000)));
	auto work = [&](int64_t a, int64_t b) {
		ptm_local_handle_t h = ptm_initialize_local();
		int8_t out_idx[PTM_MAX_INPUT_POINTS];
		double U[9], P[9], F_res[3];
		for (int64_t i = a; i < b; i++) {
			const double (*best)[3] = NULL;
			int best_i = -1;
			ptm_index(h, (size_t)i, get_neighbours, (void*)&t, flags, conventional != 0, &type[i], &alloy[i], &scale[i],
			          &rmsd[i], &q[4 * i], &F[9 * i], F_res, U, P, &inter[i], &lattice[i], &best_i, &best, out_idx);
			tmpl[i] = best_i;
		}
		ptm_uninitialize_local(h);
	};
	std::vector<std::thread> pool;
	int64_t step = (n + nthreads - 1) / nthreads;
	for (int64_t a = 0; a < n; a += step) pool.emplace_back(work, a, std::min(n, a + step));
	for (auto& th : pool) th.join();
	return 0;
}
