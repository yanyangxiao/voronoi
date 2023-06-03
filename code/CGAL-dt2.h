
#ifndef DELAUNAY_TRIANGULATION_2D_H
#define DELAUNAY_TRIANGULATION_2D_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

template <typename K, typename Vbb>
class MyDT2Vertex : public Vbb
{
public:
	typedef typename K::Point_2     Point;

	typedef typename Vbb::Triangulation_data_structure TDS;
	typedef typename TDS::Face_handle   Face_handle;
	typedef typename TDS::Vertex_handle Vertex_handle;

	template < typename TDS2 >
	struct Rebind_TDS
	{
		typedef typename Vbb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef MyDT2Vertex<K, Vb2> Other;
	};

protected:
	int _index;

public:
	MyDT2Vertex()
		: Vbb(), _index(-1)
	{
	}

	MyDT2Vertex(const Point& p)
		: Vbb(p), _index(-1)
	{
	}

	MyDT2Vertex(const Point& p, Face_handle f)
		: Vbb(p, f), _index(-1)
	{
	}

	MyDT2Vertex(const MyDT2Vertex &rhs)
		: Vbb(rhs), _index(rhs._index)
	{

	}

	~MyDT2Vertex()
	{
	}

	MyDT2Vertex& operator= (const MyDT2Vertex &rhs)
	{
		if (this != &rhs)
		{
			Vbb::operator= (rhs);
			_index = rhs._index;
		}

		return *this;
	}

	inline int index() const { return _index; }
	inline void set_index(int i) { _index = i; }
};

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel                CGALKernel;

// Vertex
typedef CGAL::Triangulation_vertex_base_2<CGALKernel>                      CGALDT2VertexBase;
typedef MyDT2Vertex<CGALKernel, CGALDT2VertexBase>                         CGALDT2Vertex;

// Face
typedef CGAL::Triangulation_face_base_2<CGALKernel>                        CGALDT2Face;

// Triangulation
typedef CGAL::Triangulation_data_structure_2<CGALDT2Vertex, CGALDT2Face>   CGALDT2DataStructure;
typedef CGAL::Delaunay_triangulation_2<CGALKernel, CGALDT2DataStructure>   CGALDelaunayTriangulation2;

class MyCGALDT2 : public CGALDelaunayTriangulation2
{
	typedef CGALDelaunayTriangulation2 baseclass;

public:
	enum {Dimension = 2};

protected:
	std::vector<Vertex_handle> _vertices;
	
public:
	inline size_t vertices_number() const
	{
		return _vertices.size(); 
	}

	inline size_t points_number() const
	{
		return _vertices.size();
	}

	void set_vertices(const double *sites, int nb)
	{
		baseclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Face_handle hf = 0;
		for (int i = 0; i < nb; ++i)
		{
			Point p(sites[2 * i], sites[2 * i + 1]);

			_vertices[i] = insert(p, hf);
			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			hf = _vertices[i]->face();			
		}
	}

	const double* position(int v) const
	{
		return &(_vertices[v]->point()[0]);
	}

	int nearest_vertex_index(const double *test) const
	{
		return nearest_vertex(Point(test[0], test[1]))->index();
	}

	int nearest_point(const double *test) const
	{
		return nearest_vertex(Point(test[0], test[1]))->index();
	}

	void vertex_neighbors(int v, std::vector<size_t> &vneighs) const
	{
		vneighs.clear();

		Vertex_circulator vvit = incident_vertices(_vertices[v]);
		Vertex_circulator vvend = vvit;
		CGAL_For_all(vvit, vvend)
		{
			int nv = vvit->index();
			if (-1 == nv)
				continue;

			vneighs.push_back(size_t(nv));
		}
	}

	template <typename Bisect>
	void vertex_bisectors(int v, std::vector<Bisect> &vbisects) const
	{
		vbisects.clear();

		const double *vp = position(v);

		Vertex_circulator vvit = incident_vertices(_vertices[v]);
		Vertex_circulator vvend = vvit;
		CGAL_For_all(vvit, vvend)
		{
			int nv = vvit->index();
			if (-1 == nv)
				continue;

			const double *q = position(nv);

			Bisect bisect;
			bisect.neigh = nv;
			bisect.plane[0] = (q[0] + vp[0]) * 0.5;
			bisect.plane[1] = (q[1] + vp[1]) * 0.5;
			bisect.plane[2] = q[0] - vp[0];
			bisect.plane[3] = q[1] - vp[1];

			vbisects.push_back(bisect);
		}
	}
};

#endif
