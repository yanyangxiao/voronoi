
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "trimesh.h"
#include "voronoi/facet.hpp"

// 1 dt method
#include "CGAL-dt3.h"
#include "voronoi/voronoi-dt.hpp"
// 2 knn method
#include "knn.hpp"
#include "voronoi/voronoi-knn.hpp"

#include "timer.h"
#include "xlog.h"

int main(int argc, char *argv[])
{
	typedef double                 MyFloat;
	typedef xyy::Facet<3, MyFloat> MyFacet3;

	typedef MyCGALDT3                                      MyDT3;
	typedef xyy::DelaunayVoronoi<MyDT3, MyFacet3, MyFloat> MyDTVoronois;

	typedef xyy::KNN<3, MyFloat>                           MyKNN3;
	typedef xyy::KnnVoronoi<MyKNN3, MyFacet3, MyFloat>     MyKnnVoronois;

	MyTrimesh tmesh;
	if (!OpenMesh::IO::read_mesh(tmesh, argv[1]))
		return 0;

	std::vector<MyFacet3> domain;
	convert_mesh<MyTrimesh, MyFacet3, MyFloat>(tmesh, domain);

	std::vector<MyFloat> sites;
	load_sites<3, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 3;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 dt voronoi:");

	MyDT3 dt;
	MyDTVoronois avoro;

	Timer t;
	t.start();
	dt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(dt, domain);
	t.stop();
	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyDTVoronois>(avoro, "voronois-dt.off");

	t.start();
	avoro.omp_compute(dt, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyDTVoronois>(avoro, "voronois-dt-omp.off");

	xlog("2 knn voronoi:");

	MyKNN3 knn;
	MyKnnVoronois bvoro;

	t.start();
	knn.set_points(&sites[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	bvoro.set_minneigh(16);

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyKnnVoronois>(bvoro, "voronois-knn.off");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyKnnVoronois>(bvoro, "voronois-knn-omp.off");

	return 0;
}
