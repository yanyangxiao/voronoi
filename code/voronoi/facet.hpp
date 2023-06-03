
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef FACET_HPP
#define FACET_HPP

#include <vector>

namespace xyy
{
    template <int DIM, typename REAL>
	class Facet
	{
	public:
		enum { Dimension = DIM };

	protected:
		std::vector<REAL> _points;
		std::vector<int>  _domneighs;
		std::vector<int>  _siteneighs;

	public:
		Facet() { }

		Facet(const Facet &rhs)
		{
			_points = rhs._points;
			_domneighs = rhs._domneighs;
			_siteneighs = rhs._siteneighs;
		}

		Facet& operator=(const Facet &rhs)
		{
			_points = rhs._points;
			_domneighs = rhs._domneighs;
			_siteneighs = rhs._siteneighs;

			return *this;
		}

		~Facet()
		{
			clear();
		}

		void clear()
		{
			_points.clear();
			_domneighs.clear();
			_siteneighs.clear();
		}

		size_t points_number() const
		{
			return _siteneighs.size();
		}

		size_t walls_number() const
		{
			return _siteneighs.size();
		}

		bool valid() const
		{
			return walls_number() > 2;
		}

		void add_point(const REAL *p, int domneigh, int siteneigh)
		{
			for (int d = 0; d < DIM; ++d)
				_points.push_back(p[d]);

			_domneighs.push_back(domneigh);
			_siteneighs.push_back(siteneigh);
		}

		REAL* point(int i)
		{
			return &_points[DIM * i];
		}

		const REAL* point(int i) const
		{
			return &_points[DIM * i];
		}

		int wall_domneigh(int i) const
		{
			return _domneighs[i];
		}

		int wall_siteneigh(int i) const
		{
			return _siteneighs[i];
		}

		template <int PlaneDIM>
		bool clip(int nv, const REAL *pp, const REAL *pn, Facet *result) const;

		void center(REAL *cent) const;
		REAL max_distance(const REAL *query, int &maxk) const;
	};

	template <int PlaneDIM, int QueryDIM, typename REAL>
	REAL plane_side(const REAL *point, const REAL *normal, const REAL *query)
	{
		REAL result = 0;
		for (int d = 0; d < QueryDIM; ++d)
			result += (query[d] - point[d]) * normal[d];

		for (int k = QueryDIM; k < PlaneDIM; ++k)
			result += (-point[k] * normal[k]);

		return result;
	}

	template <int DIM, typename REAL>
	template <int PlaneDIM>
	bool Facet<DIM, REAL>::clip(int nv, const REAL *pp, const REAL *pn, Facet *result) const
	{
		result->clear();

		bool clipped = false;

		int wallnb = (int)walls_number();
		int prev = wallnb - 1;
		const REAL *p = point(prev);
		REAL pside = plane_side<PlaneDIM, DIM, REAL>(pp, pn, p);

		for (int k = 0; k < wallnb; ++k)
		{
			const REAL *q = point(k);
			REAL qside = plane_side<PlaneDIM, DIM, REAL>(pp, pn, q);

			REAL sign = pside * qside;
			if (sign < 0.0) // intersect!
			{
				REAL x[DIM];
				REAL lambda = abs(pside) / (abs(pside) + abs(qside));
				for (int d = 0; d < DIM; ++d)
					x[d] = p[d] + lambda * (q[d] - p[d]);

				if (pside > 0.0)
				{
					result->add_point(x, wall_domneigh(prev), wall_siteneigh(prev));
				}
				else
				{
					result->add_point(p, wall_domneigh(prev), wall_siteneigh(prev));
					result->add_point(x, -1, nv);
				}

				clipped = true;
			}
			else if (sign > 0.0)
			{
				if (pside < 0.0)
				{
					result->add_point(p, wall_domneigh(prev), wall_siteneigh(prev));
				}
			}
			else // (0.0 == sign)
			{
				if (0.0 == pside)
				{
					int dn = wall_domneigh(prev);
					int sn = wall_siteneigh(prev);
					if (qside >= 0.0)
					{
						dn = -1;
						sn = nv;
					}

					result->add_point(p, dn, sn);
				}
				else // qside == 0.0
				{
					if (pside < 0.0)
					{
						result->add_point(p, wall_domneigh(prev), wall_siteneigh(prev));
					}
				}
			}

			prev = k;
			p = q;
			pside = qside;
		}

		return clipped;
	}

	template <int DIM, typename REAL>
	void Facet<DIM, REAL>::center(REAL *cent) const
	{
		int pnb = (int)points_number();
		assert(pnb > 2);

		for (int d = 0; d < DIM; ++d)
			cent[d] = REAL(0.0);

		for (int i = 0; i < pnb; ++i)
		{
			for (int d = 0; d < DIM; ++d)
				cent[d] += _points[DIM * i + d];
		}

		for (int d = 0; d < DIM; ++d)
			cent[d] /= pnb;
	}

	template <int DIM, typename REAL>
	REAL Facet<DIM, REAL>::max_distance(const REAL *query, int &maxk) const
	{
		REAL maxd = 0;
		maxk = -1;

		int pnb = (int)points_number();
		for (int k = 0; k < pnb; ++k)
		{
			const REAL *p = point(k);

			REAL dist = 0;
			for (int d = 0; d < DIM; ++d)
			{
				REAL dx = p[d] - query[d];
				dist += dx * dx;
			}
			
			if (dist > maxd)
			{
				maxd = dist;
				maxk = k;
			}
		}

		return maxd;
	}

	template <typename REAL>
	REAL area(const Facet<2, REAL> &face)
	{
		int pnb = (int)face.points_number();
		assert(pnb > 2);

		REAL area = REAL(0.0);

		const REAL *p0 = face.point(0);
		const REAL *p1 = face.point(1);

		REAL dira[2] = { p1[0] - p0[0], p1[1] - p0[1] };

		for (int i = 2; i < pnb; ++i)
		{
			const REAL *p2 = face.point(i);
			REAL dirb[2] = { p2[0] - p0[0], p2[1] - p0[1] };

			area += REAL(0.5) * (dira[0] * dirb[1] - dira[1] * dirb[0]);

			dira[0] = dirb[0];
			dira[1] = dirb[1];
		}

		return area;
	}

	template <typename REAL>
	REAL area(Facet<3, REAL> &face)
	{
		int pnb = (int)face.points_number();
		assert(pnb > 2);

		REAL area = REAL(0.0);

		const REAL *p0 = face.point(0);
		const REAL *p1 = face.point(1);

		REAL dira[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };

		for (int i = 2; i < pnb; ++i)
		{
			const REAL *p2 = face.point(i);
			REAL dirb[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };

			REAL cross[3];
			cross[0] = dira[1] * dirb[2] - dira[2] * dirb[1];
			cross[1] = dira[2] * dirb[0] - dira[0] * dirb[2];
			cross[2] = dira[0] * dirb[1] - dira[1] * dirb[0];

			area += REAL(0.5) * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

			dira[0] = dirb[0];
			dira[1] = dirb[1];
			dira[2] = dirb[2];
		}

		return area;
	}
}

#endif
