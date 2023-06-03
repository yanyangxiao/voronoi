
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef VORONOI_HPP
#define VORONOI_HPP

#include <vector>
#include <stack>
#include "bitbools.hpp"

namespace xyy
{
	template <typename Connectivity, typename SingleCell, typename REAL>
	class Voronoi
	{
	protected:
		struct StackItem
		{
			int site;      // site index
			int element;   // element index of input domain

			StackItem(int v = -1, int t = -1)
				: site(v), element(t)
			{ }
		};

		struct Bisector
		{
			int neigh;
			REAL plane[2 * Connectivity::Dimension];
			REAL radius;

			Bisector() {}
		};

	public:
		typedef std::vector<SingleCell> Cell;

	protected:
		std::vector<Cell> _cells;

	public:
		Voronoi() { }
		~Voronoi()	{ clear(); }

		void clear() { _cells.clear(); }
		size_t cells_number() const { return _cells.size(); }
		const Cell& cell(int v) const { return _cells[v]; }
		Cell& cell(int v) { return _cells[v]; }

		void compute(const Connectivity &con, const std::vector<SingleCell> &domain);

		virtual bool clip(
			const Connectivity &con,
			const SingleCell &domcell,
			int v,
			std::vector<Bisector> &vbisects,
			SingleCell &result) const = 0;

		virtual void omp_compute(const Connectivity &con, const std::vector<SingleCell> &domain) { }
	};

	template <typename Connectivity, typename SingleCell, typename REAL>
	void Voronoi<Connectivity, SingleCell, REAL>::compute(const Connectivity &con, const std::vector<SingleCell> &domain)
	{
		if (domain.empty())
			return;

		int vnb = (int)con.points_number();
		if (vnb < 3)
			return;
		
		_cells.clear();
		_cells.resize(vnb);

		std::vector<std::vector<Bisector>> bisects;
		bisects.clear();
		bisects.resize(vnb);

		int tnb = (int)domain.size();

		std::vector<BitBools> svisit(vnb, BitBools(tnb, false));
		BitBools dvisit(tnb, false);

		std::stack<StackItem> mystack;
		
		REAL query[Connectivity::Dimension] = { 0.0 };
		for (int t = 0; t < tnb; ++t)
		{
			if (dvisit.is_true(t))
				continue;

			const REAL *p = domain[t].point(0);
			
			for (int k = 0; k < SingleCell::Dimension; ++k)
				query[k] = p[k];

			int near = (int)con.nearest_point(query);

			mystack.push(StackItem(near, t));

			svisit[near].set_true(t);
			dvisit.set_true(t);

			while (!mystack.empty())
			{
				int site = mystack.top().site;
				int element = mystack.top().element;
				mystack.pop();

				SingleCell result;
				bool ok = clip(con, domain[element], site, bisects[site], result);
				if (!ok)
					continue;

				_cells[site].push_back(result);

				for (int i = 0; i < result.walls_number(); ++i)
				{
					int sn = result.wall_siteneigh(i);
					if (sn > -1 && !svisit[sn].is_true(element))
					{
						mystack.push(StackItem(sn, element));

						svisit[sn].set_true(element);
						dvisit.set_true(element);
					}

					int dn = result.wall_domneigh(i);
					if (dn > -1 && !svisit[site].is_true(dn))
					{
						mystack.push(StackItem(site, dn));

						svisit[site].set_true(dn);
						dvisit.set_true(dn);
					}
				}
			}
		}
	}
}

#endif
