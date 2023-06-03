
/**
* @author: Xiao Yanyang
* @email: yanyangxiaoxyy@gmail.com
*/

#ifndef KNN_HPP
#define KNN_HPP

#include <vector>
#include "extern/nanoflann.hpp"

namespace xyy
{
	template <int DIM, typename REAL>
	class KNN
	{
		typedef KNN<DIM, REAL> thisclass;

	public:
		enum { Dimension = DIM };

		typedef nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<REAL, thisclass>,
			thisclass, DIM> KDTree;

	protected:
		const REAL  *_points;
		size_t       _nb;
		KDTree      *_kdtree;

	public:
		KNN()
			: _points(NULL), _nb(0), _kdtree(NULL)
		{
			_kdtree = new KDTree(DIM, *this, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		}

		~KNN()
		{
			if (_kdtree)
			{
				delete _kdtree;
				_kdtree = NULL;
			}
		}

		void set_points(const REAL *points, size_t nb);
		void knn_search(const REAL *query, size_t count, std::vector<size_t> &neighs) const;
		void knn_search(const REAL *query, size_t count, std::vector<size_t> &neighs, std::vector<REAL> &dists) const;
		size_t nearest_point(const REAL *query) const;

		// access
		size_t points_number() const
		{
			return _nb;
		}

		const REAL* point(int k) const
		{
			return &_points[k * DIM];
		}

		const KDTree* kdtree() const
		{
			return _kdtree;
		}

		size_t kdtree_get_point_count() const { return _nb; }

		template <class BBOX>
		bool kdtree_get_bbox(BBOX& bb) const
		{
			return false;
		}

		REAL kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			return _points[idx * DIM + dim];
		}
	};

	template <int DIM, typename REAL>
	void KNN<DIM, REAL>::set_points(const REAL *points, size_t nb)
	{
		assert(points);

		_points = points;
		_nb = nb;

		_kdtree->buildIndex();
	}

	template <int DIM, typename REAL>
	void KNN<DIM, REAL>::knn_search(const REAL *query, size_t count, std::vector<size_t> &neighs) const
	{
		if (count > _nb)
			count = _nb;

		neighs.clear();
		neighs.resize(count);
		std::vector<REAL> dists(count);
		_kdtree->knnSearch(query, count, &neighs[0], &dists[0]);
	}

	template <int DIM, typename REAL>
	void KNN<DIM, REAL>::knn_search(const REAL *query, size_t count, std::vector<size_t> &neighs, std::vector<REAL> &dists) const
	{
		if (count > _nb)
			count = _nb;

		neighs.clear();
		neighs.resize(count);
		dists.clear();
		dists.resize(count);
		_kdtree->knnSearch(query, count, &neighs[0], &dists[0]);
	}

	template <int DIM, typename REAL>
	size_t KNN<DIM, REAL>::nearest_point(const REAL *query) const
	{
		std::vector<size_t> neighs(1);
		std::vector<REAL> dists(1);
		_kdtree->knnSearch(query, 1, &neighs[0], &dists[0]);

		return neighs[0];
	}
}

#endif
