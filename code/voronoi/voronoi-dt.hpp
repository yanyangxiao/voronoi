
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef VORONOI_DT_HPP
#define VORONOI_DT_HPP

#include "voronoi.hpp"

namespace xyy
{
	template <typename Delaunay, typename SingleCell, typename REAL>
	class DelaunayVoronoi : public Voronoi<Delaunay, SingleCell, REAL>
	{
	public:
		bool clip(
			const Delaunay &dt,
			const SingleCell &domcell,
			int v,
			std::vector<Bisector> &vbisects,
			SingleCell &result) const;
		
		void omp_compute(const Delaunay &dt, const std::vector<SingleCell> &domain);
	};

	template <typename Delaunay, typename SingleCell, typename REAL>
	bool DelaunayVoronoi<Delaunay, SingleCell, REAL>::clip(
		const Delaunay &dt,
		const SingleCell &domcell,
		int v,
		std::vector<Bisector> &vbisects,
		SingleCell &result) const
	{
		if (vbisects.empty())
			dt.vertex_bisectors<Bisector>(v, vbisects);

		if (vbisects.empty())
			return false;
		
		result = domcell;

		SingleCell tempcell;
		SingleCell *ping = &result;
		SingleCell *pong = &tempcell;

		for (int k = 0; k < (int)vbisects.size(); ++k)
		{
			ping->clip<Delaunay::Dimension>(
				vbisects[k].neigh,
				&vbisects[k].plane[0],
				&vbisects[k].plane[Delaunay::Dimension],
				pong);

			if (!pong->valid())
				return false;

			SingleCell *temp = ping;
			ping = pong;
			pong = temp;
		}

		if (ping != &result)
			result = *ping;

		return true;
	}

	template <typename Delaunay, typename SingleCell, typename REAL>
	void DelaunayVoronoi<Delaunay, SingleCell, REAL>::omp_compute(const Delaunay &dt, const std::vector<SingleCell> &domain)
	{
		if (domain.empty())
			return;

		int vnb = (int)dt.points_number();
		if (vnb < 3)
			return;

		_cells.clear();
		_cells.resize(vnb);

		std::vector<std::vector<Bisector>> bisects;
		bisects.clear();
		bisects.resize(vnb);

		int tnb = (int)domain.size();

		// getting neighbors cannot be parallel
		for (int v = 0; v < vnb; ++v)
			dt.vertex_bisectors(v, bisects[v]);

		if (tnb < 20)
		{
#pragma omp parallel for
			for (int v = 0; v < vnb; ++v)
			{
				for (int t = 0; t < tnb; ++t)
				{
					SingleCell result;
					bool ok = clip(dt, domain[t], v, bisects[v], result);
					if (!ok)
						continue;

					_cells[v].push_back(result);
				}
			}
		}
		else
		{
			// version 1
			/*std::vector<BitBools> svisit(vnb, BitBools(tnb, false));
			BitBools dvisit(tnb, false);

			std::vector<StackItem> arrA, arrB;
			std::vector<StackItem> *ping = &arrA;
			std::vector<StackItem> *pong = &arrB;

			for (int t = 0; t < tnb; ++t)
			{
				if (dvisit.is_true(t))
					continue;

				const REAL *p = domain[t].point(0);

				REAL query[Delaunay::Dimension] = { 0.0 };
				for (int k = 0; k < SingleCell::Dimension; ++k)
					query[k] = p[k];

				int near = (int)dt.nearest_point(query);

				ping->push_back(StackItem(near, t));

				svisit[near].set_true(t);
				dvisit.set_true(t);

				while (!ping->empty())
				{
					pong->clear();

					int nb = (int)ping->size();
#pragma omp parallel for
					for (int k = 0; k < nb; ++k)
					{
						int site = (*ping)[k].site;
						int element = (*ping)[k].element;

						SingleCell result;
						bool ok = clip(dt, domain[element], site, bisects[site], result);
						if (!ok)
							continue;
						
#pragma omp critical
						{
							_cells[site].push_back(result);
						}

						for (int i = 0; i < result.walls_number(); ++i)
						{
							int sn = result.wall_siteneigh(i);
							int dn = result.wall_domneigh(i);
#pragma omp critical
							{
								if (sn > -1 && !svisit[sn].is_true(element))
								{
									pong->push_back(StackItem(sn, element));

									svisit[sn].set_true(element);
									dvisit.set_true(element);
								}
								
								if (dn > -1 && !svisit[site].is_true(dn))
								{
									pong->push_back(StackItem(site, dn));

									svisit[site].set_true(dn);
									dvisit.set_true(dn);
								}
							}
						}
					}

					std::vector<StackItem> *temp = ping;
					ping = pong;
					pong = temp;
				}
			}*/

			// version 2
			std::vector<std::vector<bool>> visited(vnb, std::vector<bool>(tnb, false));
			std::vector<std::stack<int>> stacks(vnb);

#pragma omp parallel for
			for (int t = 0; t < tnb; ++t)
			{
				const REAL *p = domain[t].point(0);

				REAL query[Delaunay::Dimension] = { 0.0 };
				for (int k = 0; k < SingleCell::Dimension; ++k)
					query[k] = p[k];

				int near = (int)dt.nearest_point(query); // can be parallel?

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
					bool ok = clip(dt, domain[element], v, bisects[v], result);
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
