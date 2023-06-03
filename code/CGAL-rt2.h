
#ifndef REGULAR_TRIANGULATION_2D_H
#define REGULAR_TRIANGULATION_2D_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

template <typename K, typename Vbb>
class MyRT2Vertex : public Vbb
{
public:
	typedef typename Vbb::Bare_point                   Bare_point;
	typedef typename Vbb::Weighted_point               Weighted_point;

	typedef typename Vbb::Triangulation_data_structure TDS;
	typedef typename TDS::Face_handle                  Face_handle;
	typedef typename TDS::Vertex_handle                Vertex_handle;

	template < typename TDS2 >
	struct Rebind_TDS
	{
		typedef typename Vbb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef MyRT2Vertex<K, Vb2> Other;
	};

protected:
	int _index;

public:
	MyRT2Vertex()
		: Vbb(), _index(-1)
	{
	}

	MyRT2Vertex(const Weighted_point& p)
		: Vbb(p), _index(-1)
	{
	}

	MyRT2Vertex(const Weighted_point& p, Face_handle f)
		: Vbb(p, f), _index(-1)
	{
	}

	MyRT2Vertex(const MyRT2Vertex &rhs)
		: Vbb(rhs), _index(rhs._index)
	{

	}

	~MyRT2Vertex()
	{
	}

	MyRT2Vertex& operator= (const MyRT2Vertex &rhs)
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
typedef CGAL::Exact_predicates_inexact_constructions_kernel                Kernel;

// Vertex
typedef CGAL::Regular_triangulation_vertex_base_2<Kernel>                  CGALRT2VertexBase;
typedef MyRT2Vertex<Kernel, CGALRT2VertexBase>                             CGALRT2Vertex;

// Face
typedef CGAL::Regular_triangulation_face_base_2<Kernel>                    CGALRT2Face;

// Triangulation
typedef CGAL::Triangulation_data_structure_2<CGALRT2Vertex, CGALRT2Face>   CGALRT2DataStructure;
typedef CGAL::Regular_triangulation_2<Kernel, CGALRT2DataStructure>        CGALRegularTriangulation2;


class MyCGALRT2 : public CGALRegularTriangulation2
{
	typedef CGALRegularTriangulation2 baseclass;

public:
	enum { Dimension = 2 };

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

	void set_vertices(const double *sites, const double *weights, int nb)
	{
		baseclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Face_handle hf = 0;
		for (int i = 0; i < nb; ++i)
		{
			Weighted_point wp(Bare_point(sites[2 * i], sites[2 * i + 1]), weights[i]);

			_vertices[i] = insert(wp, hf);
			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			hf = _vertices[i]->face();
		}
	}

	void set_vertices(const double *pw, int nb)
	{
		baseclass::clear();
		_vertices.clear();
		_vertices.resize(nb, 0);

		Face_handle hf = 0;
		for (int i = 0; i < nb; ++i)
		{
			Weighted_point wp(Bare_point(pw[3 * i], pw[3 * i + 1]), pw[3 * i + 2]);

			_vertices[i] = insert(wp, hf);
			if (_vertices[i] == 0)
				continue;

			_vertices[i]->set_index(i);

			hf = _vertices[i]->face();
		}
	}

	const double* position(int v) const
	{
		return &(_vertices[v]->point().point()[0]);
	}

	double weight(int v) const
	{
		return _vertices[v]->point().weight();
	}

	int nearest_vertex_index(const double *test) const
	{
		return nearest_power_vertex(Bare_point(test[0], test[1]))->index();
	}

	int nearest_point(const double *test) const
	{
		return nearest_power_vertex(Bare_point(test[0], test[1]))->index();
	}

	void vertex_neighbors(int v, std::vector<size_t> &vneighs) const
	{
		vneighs.clear();

		if (_vertices[v] == 0 || v != _vertices[v]->index() || _vertices[v]->is_hidden())
			return;

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

	void middle_point(const double *p, double pw, const double *q, double qw, double *mid) const
	{
		double dx = q[0] - p[0];
		double dy = q[1] - p[1];
		double len = std::sqrt(dx * dx + dy * dy);

		dx /= len;
		dy /= len;

		double dist = 0.5 * (len + (pw - qw) / len);

		mid[0] = p[0] + dist * dx;
		mid[1] = p[1] + dist * dy;
	}

	template <typename Bisect>
	void vertex_bisectors(int v, std::vector<Bisect> &vbisects) const
	{
		vbisects.clear();

		if (_vertices[v] == 0 || v != _vertices[v]->index() || _vertices[v]->is_hidden())
			return;

		const double *vp = position(v);
		double vw = weight(v);

		Vertex_circulator vvit = incident_vertices(_vertices[v]);
		Vertex_circulator vvend = vvit;
		CGAL_For_all(vvit, vvend)
		{
			int nv = vvit->index();
			if (-1 == nv)
				continue;

			const double *q = position(nv);
			double w = weight(nv);

			double mid[2];
			middle_point(vp, vw, q, w, mid);

			Bisect bisect;
			bisect.neigh = nv;
			bisect.plane[0] = mid[0];
			bisect.plane[1] = mid[1];
			bisect.plane[2] = q[0] - vp[0];
			bisect.plane[3] = q[1] - vp[1];

			vbisects.push_back(bisect);
		}
	}
};

#endif
