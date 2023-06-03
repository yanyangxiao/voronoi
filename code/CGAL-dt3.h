
#ifndef DELAUNAY_TRIANGULATION_3D_H
#define DELAUNAY_TRIANGULATION_3D_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
class MyDT3Vertex : public Vb
{
	typedef Vb baseclass;

public:
	typedef typename Vb::Cell_handle        Cell_handle;
	typedef typename Vb::Point              Point;

	template < typename TDS2 >
	struct Rebind_TDS
	{
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef MyDT3Vertex<Gt, Vb2>   Other;
	};

protected:
	int _index;

public:
	MyDT3Vertex()
		: baseclass(), _index(-1)
	{}

	MyDT3Vertex(const Point &p)
		: baseclass(p), _index(-1)
	{}

	MyDT3Vertex(const Point &p, Cell_handle c)
		: baseclass(p, c), _index(-1)
	{}

	MyDT3Vertex(Cell_handle c)
		: baseclass(c), _index(-1)
	{}

	MyDT3Vertex(const MyDT3Vertex &rhs)
		: baseclass(rhs), _index(rhs._index)
	{}

	MyDT3Vertex& operator= (const MyDT3Vertex &rhs)
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
typedef CGAL::Triangulation_vertex_base_3<CGALKernel>                    CGALDT3VertexBase;
typedef MyDT3Vertex<CGALKernel, CGALDT3VertexBase>                       CGALDT3Vertex;

typedef CGAL::Delaunay_triangulation_cell_base_3<CGALKernel>             CGALDT3Cell;

typedef CGAL::Triangulation_data_structure_3<CGALDT3Vertex, CGALDT3Cell> CGALDT3DataStructure;
typedef CGAL::Delaunay_triangulation_3<CGALKernel, CGALDT3DataStructure> CGALDelaunayTriangulation3;

class MyCGALDT3 : public CGALDelaunayTriangulation3
{
	typedef CGALDelaunayTriangulation3 superclass;

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
		superclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Cell_handle ch;
		for (int i = 0; i < nb; ++i)
		{
			Point p(sites[3 * i], sites[3 * i + 1], sites[3 * i + 2]);
			_vertices[i] = insert(p, ch);

			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			ch = _vertices[i]->cell();
		}
	}

	inline int nearest_vertex_index(const double *test) const
	{
		return nearest_vertex(Point(test[0], test[1], test[2]))->index();
	}

	inline int nearest_point(const double *test) const
	{
		return nearest_vertex(Point(test[0], test[1], test[2]))->index();
	}

	const double* position(int v) const
	{
		return &(_vertices[v]->point()[0]);
	}

	void vertex_neighbors(int v, std::vector<size_t> &vneighs) const
	{
		vneighs.clear();

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

	template <typename Bisect>
	void vertex_bisectors(int v, std::vector<Bisect> &vbisects) const
	{
		vbisects.clear();

		std::vector<Vertex_handle> adjacents;
		finite_adjacent_vertices(_vertices[v], std::back_inserter(adjacents));

		const double *vp = position(v);

		for (int k = 0; k < adjacents.size(); ++k)
		{
			int nv = adjacents[k]->index();
			if (-1 == nv)
				continue;

			const double *q = position(nv);

			Bisect bisect;
			bisect.neigh = nv;
			bisect.plane[0] = (q[0] + vp[0]) * 0.5;
			bisect.plane[1] = (q[1] + vp[1]) * 0.5;
			bisect.plane[2] = (q[2] + vp[2]) * 0.5;
			bisect.plane[3] = q[0] - vp[0];
			bisect.plane[4] = q[1] - vp[1];
			bisect.plane[5] = q[2] - vp[2];

			vbisects.push_back(bisect);
		}
	}
};

#endif
