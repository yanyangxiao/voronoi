
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef POLYMESH_HPP
#define POLYMESH_HPP

#include <vector>

namespace xyy
{
	template <typename REAL>
	class Polymesh
	{
	public:
		enum { Dimension = 3 };

	protected:
		std::vector<int>   _facets;
		std::vector<int>   _domneighs;
		std::vector<int>   _siteneighs;
		std::vector<REAL>  _points;

	public:
		Polymesh()
		{
			_facets.push_back(0);
		}

		Polymesh(const Polymesh &rhs)
		{
			_facets = rhs._facets;
			_domneighs = rhs._domneighs;
			_siteneighs = rhs._siteneighs;
			_points = rhs._points;
		}

		Polymesh& operator=(const Polymesh &rhs)
		{
			_facets = rhs._facets;
			_domneighs = rhs._domneighs;
			_siteneighs = rhs._siteneighs;
			_points = rhs._points;

			return *this;
		}

		~Polymesh()
		{
			clear();
		}

		void clear()
		{
			_facets.clear();
			_domneighs.clear();
			_siteneighs.clear();
			_points.clear();

			_facets.push_back(0);
		}

		size_t facets_number() const
		{
			return _siteneighs.size();
		}

		size_t walls_number() const
		{
			return _siteneighs.size();
		}

		size_t points_number() const
		{
			return _points.size() / 3;
		}

		bool valid() const
		{
			return walls_number() > 2;
		}

		void begin_facet()
		{ }

		void add_point(const REAL *p)
		{
			_points.push_back(p[0]);
			_points.push_back(p[1]);
			_points.push_back(p[2]);
		}

		void add_point(REAL x, REAL y, REAL z)
		{
			_points.push_back(x);
			_points.push_back(y);
			_points.push_back(z);
		}

		void end_facet(int domneigh, int siteneigh)
		{
			_facets.push_back((int)points_number());
			_domneighs.push_back(domneigh);
			_siteneighs.push_back(siteneigh);
		}

		int facet_begin(int f) const
		{
			return _facets[f];
		}

		int facet_end(int f) const
		{
			return _facets[f + 1];
		}

		int facet_size(int f) const
		{
			return facet_end(f) - facet_begin(f);
		}

		REAL* point(int i)
		{
			return &_points[3 * i];
		}

		const REAL* point(int i) const
		{
			return &_points[3 * i];
		}
		
		int wall_domneigh(int f) const
		{
			return _domneighs[f];
		}

		int wall_siteneigh(int f) const
		{
			return _siteneighs[f];
		}

		int next_around_facet(int f, int i) const
		{
			return (i + 1 == facet_end(f) ? facet_begin(f) : i + 1);
		}

		int prev_around_facet(int f, int i) const
		{
			return (i == facet_begin(f) ? facet_end(f) - 1 : i - 1);
		}

		template <int PlaneDIM>
		bool clip(int nv, const REAL *pp, const REAL *pn, Polymesh *result) const;
		
		void center(REAL *cent) const;
		void facet_center(int f, REAL *fcent) const;
		REAL facet_area(int f) const;

		REAL max_distance(const REAL *query, int &maxk) const;
	};

	template <int PlaneDIM, typename REAL>
	REAL plane_side3(const REAL *point, const REAL *normal, const REAL *query)
	{
		REAL result = (query[0] - point[0]) * normal[0] +
			(query[1] - point[1]) * normal[1] +
			(query[2] - point[2]) * normal[2];
		
		for (int k = 3; k < PlaneDIM; ++k)
			result += (-point[k] * normal[k]);

		return result;
	}

	template <typename REAL>
	template <int PlaneDIM>
	bool Polymesh<REAL>::clip(int nv, const REAL *pp, const REAL *pn, Polymesh *result) const
	{
		result->clear();

		bool clipped = false;

		int wallnb = (int)walls_number();
		
		std::vector<int> intersects(2 * wallnb, -1);

		for (int f = 0; f < wallnb; ++f)
		{
			int first = facet_begin(f);
			int prev = prev_around_facet(f, first);
			const REAL *p = point(prev);
			REAL pside = plane_side3<PlaneDIM>(pp, pn, p);

			bool hasbegan = false;
			for (int fv = first; fv < facet_end(f); ++fv)
			{
				const REAL *q = point(fv);
				REAL qside = plane_side3<PlaneDIM>(pp, pn, q);

				REAL sign = pside * qside;
				if (sign < 0.0) // intersect!
				{
					REAL x[3];
					REAL lambda = abs(pside) / (abs(pside) + abs(qside));
					x[0] = p[0] + lambda * (q[0] - p[0]);
					x[1] = p[1] + lambda * (q[1] - p[1]);
					x[2] = p[2] + lambda * (q[2] - p[2]);

					if (!hasbegan)
					{
						hasbegan = true;
						result->begin_facet();
					}

					if (pside > 0.0)
					{
						intersects[2 * f] = (int)result->points_number();
						result->add_point(x);
					}
					else
					{
						result->add_point(p);
						intersects[2 * f + 1] = (int)result->points_number();
						result->add_point(x);
					}

					clipped = true;
				}
				else if (sign > 0.0)
				{
					if (pside < 0.0)
					{
						if (!hasbegan)
						{
							hasbegan = true;
							result->begin_facet();
						}

						result->add_point(p);
					}
				}
				else // (0.0 == sign)
				{
					if (0.0 == pside)
					{
						if (!hasbegan)
						{
							hasbegan = true;
							result->begin_facet();
						}

						result->add_point(p);
					}
					else // qside == 0.0
					{
						if (pside < 0.0)
						{
							if (!hasbegan)
							{
								hasbegan = true;
								result->begin_facet();
							}

							result->add_point(p);
						}
					}
				}

				prev = fv;
				p = q;
				pside = qside;
			}

			if (hasbegan)
			{
				result->end_facet(wall_domneigh(f), wall_siteneigh(f));
			}
		}

		if (!clipped)
			return false;

		// connect intersections

		int idx = -1;
		for (int f = 0; f < wallnb; ++f)
		{
			if (intersects[2 * f] > -1)
			{
				idx = f;
				break;
			}
		}
		if (idx < 0)
			return false;

		std::vector<bool> visited(wallnb, false);

		result->begin_facet();
		while (idx != -1)
		{
			const REAL *p = result->point(intersects[2 * idx]);
			// fix a bug here
			// result->add_point(p);
			REAL x = p[0];
			REAL y = p[1];
			REAL z = p[2];
			result->add_point(x, y, z);
			visited[idx] = true;

			p = result->point(intersects[2 * idx + 1]);

			int next = -1;
			REAL minDist = 1e10;

			for (int f = 0; f < wallnb; ++f)
			{
				if (f == idx)
					continue;

				if (-1 == intersects[2 * f])
					continue;

				if (visited[f])
					continue;

				const REAL *q = result->point(intersects[2 * f]);

				REAL dx = q[0] - p[0];
				REAL dy = q[1] - p[1];
				REAL dz = q[2] - p[2];

				REAL dist = dx * dx + dy * dy + dz * dz;
				if (dist < minDist)
				{
					minDist = dist;
					next = f;
				}
			}
			
			idx = next;
		}
		result->end_facet(-1, nv);

		return true;
	}

	template <typename REAL>
	REAL Polymesh<REAL>::facet_area(int f) const
	{
		REAL area = REAL(0.0);

		int start = facet_begin(f);
		const REAL *p0 = point(start);

		int next = next_around_facet(f, start);
		int nextnext = next_around_facet(f, next);

		const REAL *p1 = point(next);
		REAL dira[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };

		for (int i = nextnext; i != start; i = next_around_facet(f, i))
		{
			const REAL *p2 = point(i);
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

	template <typename REAL>
	void Polymesh<REAL>::facet_center(int f, REAL *fcent) const
	{
		int pnb = facet_size(f);

		fcent[0] = REAL(0.0);
		fcent[1] = REAL(0.0);
		fcent[2] = REAL(0.0);

		for (int i = facet_begin(f); i < facet_end(f); ++i)
		{
			fcent[0] += _points[3 * i];
			fcent[1] += _points[3 * i + 1];
			fcent[2] += _points[3 * i + 2];
		}

		fcent[0] /= pnb;
		fcent[1] /= pnb;
		fcent[2] /= pnb;
	}

	template <typename REAL>
	void Polymesh<REAL>::center(REAL *cent) const
	{
		cent[0] = REAL(0.0);
		cent[1] = REAL(0.0);
		cent[2] = REAL(0.0);

		REAL totalarea = 0;

		int fnb = (int)facets_number();
		for (int f = 0; f < fnb; ++f)
		{
			REAL fcent[3];
			facet_center(f, fcent);

			REAL farea = facet_area(f);

			cent[0] += fcent[0] * farea;
			cent[1] += fcent[1] * farea;
			cent[2] += fcent[2] * farea;

			totalarea += farea;
		}

		cent[0] /= totalarea;
		cent[1] /= totalarea;
		cent[2] /= totalarea;
	}

	template <typename REAL>
	REAL Polymesh<REAL>::max_distance(const REAL *query, int &maxk) const
	{
		REAL maxd = 0;
		maxk = -1;

		int pnb = (int)points_number();
		for (int k = 0; k < pnb; ++k)
		{
			const REAL *p = point(k);

			REAL dist = 0;
			for (int d = 0; d < 3; ++d)
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
}

#endif
