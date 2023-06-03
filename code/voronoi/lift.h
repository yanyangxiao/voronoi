
#ifndef LIFT_H
#define LIFT_H

#include <vector>

template <int DIM, typename REAL>
void lift_sites(const std::vector<REAL> &sites, std::vector<REAL> &lifts, REAL times = 1)
{
	lifts.clear();
	lifts.resize(sites.size());

	memcpy(&lifts[0], &sites[0], sizeof(REAL) * sites.size());

	int vnb = (int)sites.size() / DIM;
	REAL wmax = sites[DIM - 1];
	for (int v = 1; v < vnb; ++v)
	{
		if (sites[DIM * v + DIM - 1] > wmax)
			wmax = sites[DIM * v + DIM - 1];
	}

	wmax += 1e-12;
	REAL eta = wmax * times;

	for (int v = 0; v < vnb; ++v)
	{
		lifts[DIM * v + DIM - 1] = std::sqrt(eta - sites[DIM * v + DIM - 1]);
	}
}

template <int DIM, typename REAL>
void lift_sites(const std::vector<REAL> &sites, const std::vector<REAL> &weights, std::vector<REAL> &lifts, REAL times = 1)
{
	lifts.clear();
	lifts.resize(sites.size() + weights.size());

	int vnb = (int)weights.size();
	REAL wmax = weights[0];
	for (int v = 1; v < vnb; ++v)
	{
		if (weights[v] > wmax)
			wmax = weights[v];
	}

	wmax += 1e-12;
	REAL eta = wmax * times;

	for (int v = 0; v < vnb; ++v)
	{
		for (int d = 0; d < DIM - 1; ++d)
			lifts[DIM * v + d] = sites[(DIM - 1) * v + d];

		lifts[DIM * v + DIM - 1] = std::sqrt(eta - weights[v]);
	}
}

#endif
