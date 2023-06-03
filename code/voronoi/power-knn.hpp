
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef POWER_KNN_HPP
#define POWER_KNN_HPP

#include "voronoi.hpp"

namespace xyy
{
	template <typename KNN, typename SingleCell, typename REAL>
	class KnnPower : public Voronoi<KNN, SingleCell, REAL>
	{
		typedef Voronoi<KNN, SingleCell, REAL> baseclass;

	protected:
		int _minneigh;

	public:
		KnnPower() : baseclass(), _minneigh(16) {}

		void set_minneigh(int n) { _minneigh = n; }

		bool clip(
			const KNN &knn,
			const SingleCell &domcell,
			int v,
			std::vector<Bisector> &vbisects,
			SingleCell &result) const;

		void omp_compute(const KNN &knn, const std::vector<SingleCell> &domain);
	};

	template <typename KNN, typename SingleCell, typename REAL>
	bool KnnPower<KNN, SingleCell, REAL>::clip(
		const KNN &knn,
		const SingleCell &domcell, 
		int v,
		std::vector<Bisector> &vbisects,
		SingleCell &result) const
	{
		result.clear();
		
		const REAL *vp = knn.point(v);
		REAL twin[KNN::Dimension];
		for (int d = 0; d < KNN::Dimension - 1; ++d)
			twin[d] = vp[d];

		twin[KNN::Dimension - 1] = -vp[KNN::Dimension - 1];

		int count = (int)vbisects.size();
		if (count < _minneigh)
		{
			count = _minneigh;
			std::vector<size_t> vneigh;
			knn.knn_search(twin, count, vneigh);

			vbisects.clear();
			vbisects.resize(vneigh.size());

			for (int i = 0; i < (int)vneigh.size(); ++i)
			{
				int nv = (int)vneigh[i];
				const REAL *nvpos = knn.point(nv);

				vbisects[i].neigh = nv;
				vbisects[i].radius = 0;
				for (int d = 0; d < KNN::Dimension; ++d)
				{
					vbisects[i].plane[d] = (nvpos[d] + twin[d]) * 0.5;
					REAL dx = nvpos[d] - twin[d];
					vbisects[i].plane[KNN::Dimension + d] = dx;
					vbisects[i].radius += dx * dx;
				}
				vbisects[i].radius /= 4;
			}
		}

		if (vbisects.empty())
			return false;

		result = domcell;

		SingleCell tempcell;
		SingleCell *ping = &result;
		SingleCell *pong = &tempcell;

		int maxk = -1;
		REAL maxDist = ping->max_distance(twin, maxk);
		for (int d = SingleCell::Dimension; d < KNN::Dimension; ++d)
			maxDist += twin[d] * twin[d];

		int k = 0;
		while (k < knn.points_number())
		{
			bool allin = false;
			for (; k < (int)vbisects.size(); ++k)
			{
				if (vbisects[k].neigh == v)
					continue;

				if (vbisects[k].radius > maxDist)
				{
					allin = true;
					break;
				}

				bool clipped = ping->clip<KNN::Dimension>(
					vbisects[k].neigh,
					&vbisects[k].plane[0],
					&vbisects[k].plane[KNN::Dimension],
					pong);

				if (!pong->valid())
					return false;

				SingleCell *temp = ping;
				ping = pong;
				pong = temp;

				if (clipped)
				{
					maxDist = ping->max_distance(twin, maxk);
					for (int d = SingleCell::Dimension; d < KNN::Dimension; ++d)
						maxDist += twin[d] * twin[d];
				}
			}

			if (allin || knn.points_number() <= (int)vbisects.size())
				break;

			count *= 2;
			std::vector<size_t> vneigh;
			knn.knn_search(twin, count, vneigh);

			for (int i = k; i < (int)vneigh.size(); ++i)
			{
				int nv = (int)vneigh[i];
				const REAL *nvpos = knn.point(nv);

				Bisector bisect;
				bisect.neigh = nv;
				bisect.radius = 0;
				for (int d = 0; d < KNN::Dimension; ++d)
				{
					bisect.plane[d] = (nvpos[d] + twin[d]) * 0.5;
					REAL dx = nvpos[d] - twin[d];
					bisect.plane[KNN::Dimension + d] = dx;
					bisect.radius += dx * dx;
				}
				bisect.radius /= 4;

				vbisects.push_back(bisect);
			}
		}

		if (ping != &result)
			result = *ping;

		return true;
	}

	template <typename KNN, typename SingleCell, typename REAL>
	void KnnPower<KNN, SingleCell, REAL>::omp_compute(const KNN &knn, const std::vector<SingleCell> &domain)
	{
		if (domain.empty())
			return;

		int vnb = (int)knn.points_number();
		if (vnb < 2)
			return;

		_cells.clear();
		_cells.resize(vnb);

		std::vector<std::vector<Bisector>> bisects;
		bisects.clear();
		bisects.resize(vnb);

#pragma omp parallel for
		for (int v = 0; v < vnb; ++v)
		{
			const REAL *vp = knn.point(v);
			REAL twin[KNN::Dimension];
			for (int d = 0; d < KNN::Dimension - 1; ++d)
				twin[d] = vp[d];

			twin[KNN::Dimension - 1] = -vp[KNN::Dimension - 1];

			std::vector<size_t> vneigh;
			knn.knn_search(twin, _minneigh, vneigh);

			bisects[v].clear();
			bisects[v].resize(vneigh.size());

			for (int i = 0; i < (int)vneigh.size(); ++i)
			{
				int nv = (int)vneigh[i];
				const REAL *nvpos = knn.point(nv);

				bisects[v][i].neigh = nv;
				bisects[v][i].radius = 0;
				for (int d = 0; d < KNN::Dimension; ++d)
				{
					bisects[v][i].plane[d] = (nvpos[d] + twin[d]) * 0.5;
					REAL dx = nvpos[d] - twin[d];
					bisects[v][i].plane[KNN::Dimension + d] = dx;
					bisects[v][i].radius += dx * dx;
				}
				bisects[v][i].radius /= 4;
			}
		}

		int tnb = (int)domain.size();

		if (tnb < 20)
		{
#pragma omp parallel for
			for (int v = 0; v < vnb; ++v)
			{
				for (int t = 0; t < tnb; ++t)
				{
					SingleCell result;
					bool ok = clip(knn, domain[t], v, bisects[v], result);
					if (!ok)
						continue;

					_cells[v].push_back(result);
				}
			}
		}
		else
		{
			std::vector<std::vector<bool>> visited(vnb, std::vector<bool>(tnb, false));
			std::vector<std::stack<int>> stacks(vnb);

#pragma omp parallel for
			for (int t = 0; t < tnb; ++t)
			{
				const REAL *p = domain[t].point(0);

				REAL query[KNN::Dimension] = { 0.0 };
				for (int k = 0; k < SingleCell::Dimension; ++k)
					query[k] = p[k];

				int near = (int)knn.nearest_point(query);
				
#pragma omp critical
				{
					stacks[near].push(t);
				}

				visited[near][t] = true;
			}

			bool done = false;
			while (!done)
			{
				done = true;

#pragma omp parallel for
				for (int v = 0; v < vnb; ++v)
				{
					if (stacks[v].empty())
						continue;

					int element = stacks[v].top();
					stacks[v].pop();

					SingleCell result;
					bool ok = clip(knn, domain[element], v, bisects[v], result);
					if (!ok)
						continue;

					_cells[v].push_back(result);

					for (int i = 0; i < result.walls_number(); ++i)
					{
						int sn = result.wall_siteneigh(i);
						int dn = result.wall_domneigh(i);
#pragma omp critical
						{
							if (sn > -1 && !visited[sn][element])
							{
								stacks[sn].push(element);
								visited[sn][element] = true;
							}

							if (dn > -1 && !visited[v][dn])
							{
								stacks[v].push(dn);
								visited[v][dn] = true;
							}
						}
					}

					if (!stacks[v].empty())
						done = false;
				}
			}
		}
	}	
}

#endif
