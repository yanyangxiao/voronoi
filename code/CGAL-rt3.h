
#ifndef REGULAR_TRIANGULATION_3D_H
#define REGULAR_TRIANGULATION_3D_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_3.h>

template <class Gt, class Vb = CGAL::Regular_triangulation_vertex_base_3<Gt> >
class MyRT3Vertex : public Vb
{
	typedef Vb baseclass;

public:
	typedef typename Vb::Cell_handle        Cell_handle;
	typedef typename Gt::Point_3            Bare_point;
	typedef typename Gt::Weighted_point_3   Weighted_point;

	template < typename TDS2 >
	struct Rebind_TDS
	{
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef MyRT3Vertex<Gt, Vb2> Other;
	};

protected:
	int _index;

public:
	MyRT3Vertex()
		: baseclass(), _index(-1)
	{}

	MyRT3Vertex(const Weighted_point &p)
		: baseclass(p), _index(-1)
	{}

	MyRT3Vertex(const Weighted_point &p, Cell_handle c)
		: baseclass(p, c), _index(-1)
	{}

	MyRT3Vertex(Cell_handle c)
		: baseclass(c), _index(-1)
	{}

	MyRT3Vertex(const MyRT3Vertex &rhs)
		: baseclass(rhs), _index(rhs._index)
	{}

	MyRT3Vertex& operator= (const MyRT3Vertex &rhs)
	{
		baseclass::operator= (rhs);
		_index = rhs._index;

		return *this;
	}

	inline int index() const
	{
		return _index;
	}

	inline void set_index(int i)
	{
		_index = i;
	}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel              CGALKernel;
typedef CGAL::Regular_triangulation_vertex_base_3<CGALKernel>            CGALRT3VertexBase;
typedef MyRT3Vertex<CGALKernel, CGALRT3VertexBase>                       CGALRT3Vertex;

typedef CGAL::Regular_triangulation_cell_base_3<CGALKernel>              CGALRT3Cell;

typedef CGAL::Triangulation_data_structure_3<CGALRT3Vertex, CGALRT3Cell> CGALRT3DataStructure;
typedef CGAL::Regular_triangulation_3<CGALKernel, CGALRT3DataStructure>  CGALRegularTriangulation3;

class MyCGALRT3 : public CGALRegularTriangulation3
{
	typedef CGALRegularTriangulation3 baseclass;

public:
	enum { Dimension = 3 };

protected:
	std::vector<Vertex_handle> _vertices;

public:
	inline int vertices_number() const
	{
		return (int)_vertices.size();
	}

	inline int points_number() const
	{
		return (int)_vertices.size();
	}

	void set_vertices(const double *sites, int nb)
	{
		baseclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Cell_handle ch;
		for (int i = 0; i < nb; ++i)
		{
			Weighted_point p(Bare_point(sites[4 * i], sites[4 * i + 1], sites[4 * i + 2]), sites[4 * i + 3]);
			_vertices[i] = insert(p, ch);

			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			ch = _vertices[i]->cell();
		}
	}

	void set_vertices(const double *sites, const double *weights, int nb)
	{
		baseclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Cell_handle ch;
		for (int i = 0; i < nb; ++i)
		{
			Weighted_point p(Bare_point(sites[3 * i], sites[3 * i + 1], sites[3 * i + 2]), weights[i]);
			_vertices[i] = insert(p, ch);

			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			ch = _vertices[i]->cell();
		}
	}

	inline int nearest_vertex_index(const double *test) const
	{
		return nearest_power_vertex(Bare_point(test[0], test[1], test[2]))->index();
	}

	inline int nearest_point(const double *test) const
	{
		return nearest_power_vertex(Bare_point(test[0], test[1], test[2]))->index();
	}

	const double* position(int v) const
	{
		return &(_vertices[v]->point()[0]);
	}

	double weight(int v) const
	{
		return _vertices[v]->point().weight();
	}

	void vertex_neighbors(int v, std::vector<int> &vneighs) const
	{
		vneighs.clear();
		
		if (_vertices[v] == 0 || v != _vertices[v]->index())
			return;

		const Cell_handle& c = _vertices[v]->cell();
		if (!c->has_vertex(_vertices[v]))
			return;

		std::vector<Vertex_handle> adjacents;
		finite_adjacent_vertices(_vertices[v], std::back_inserter(adjacents));

		for (int k = 0; k < adjacents.size(); ++k)
		{
			int nv = adjacents[k]->index();
			if (-1 == nv)
				continue;

			vneighs.push_back(nv);
		}
	}

	void middle_point(const double *p, double pw, const double *q, double qw, double *mid) const
	{
		double dx = q[0] - p[0];
		double dy = q[1] - p[1];
		double dz = q[2] - p[2];
		double len = std::sqrt(dx * dx + dy * dy + dz * dz);

		dx /= len;
		dy /= len;
		dz /= len;

		double dist = 0.5 * (len + (pw - qw) / len);

		mid[0] = p[0] + dist * dx;
		mid[1] = p[1] + dist * dy;
		mid[2] = p[2] + dist * dz;
	}

	template <typename Bisect>
	void vertex_bisectors(int v, std::vector<Bisect> &vbisects) const
	{
		vbisects.clear();

		if (_vertices[v] == 0 || v != _vertices[v]->index())
			return;

		const Cell_handle& c = _vertices[v]->cell();
		if (!c->has_vertex(_vertices[v]))
			return;

		const double *vp = position(v);
		double vw = weight(v);

		std::vector<Vertex_handle> adjacents;
		finite_adjacent_vertices(_vertices[v], std::back_inserter(adjacents));

		for (int k = 0; k < adjacents.size(); ++k)
		{
			int nv = adjacents[k]->index();
			if (-1 == nv)
				continue;

			const double *q = position(nv);
			double w = weight(nv);

			double mid[3];
			middle_point(vp, vw, q, w, mid);

			Bisect bisect;
			bisect.neigh = nv;
			bisect.plane[0] = mid[0];
			bisect.plane[1] = mid[1];
			bisect.plane[2] = mid[2];
			bisect.plane[3] = q[0] - vp[0];
			bisect.plane[4] = q[1] - vp[1];
			bisect.plane[5] = q[2] - vp[2];

			vbisects.push_back(bisect);
		}
	}
};

#endif
