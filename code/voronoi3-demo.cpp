
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "voronoi/polymesh.hpp"

// dt
#include "CGAL-dt3.h"
#include "voronoi/voronoi-dt.hpp"
// knn
#include "knn.hpp"
#include "voronoi/voronoi-knn.hpp"

#include "timer.h"
#include "xlog.h"

int main(int argc, char *argv[])
{
	typedef double                    MyFloat;
	typedef xyy::Polymesh<MyFloat>    MyPolymesh;

	typedef MyCGALDT3                                        MyDT3;
	typedef xyy::DelaunayVoronoi<MyDT3, MyPolymesh, MyFloat> MyDTVoronoi3;

	typedef xyy::KNN<3, MyFloat>                             MyKNN3;
	typedef xyy::KnnVoronoi<MyKNN3, MyPolymesh, MyFloat>     MyKnnVoronoi3;

	std::vector<MyPolymesh> domain;
	load_tets<MyPolymesh, MyFloat>(argv[1], domain);

	std::vector<MyFloat> sites;
	load_sites<3, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 3;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 dt voronoi:");

	MyDT3 dt;
	MyDTVoronoi3 avoro;

	Timer t;
	t.start();
	dt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(dt, domain);
	t.stop();
	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyDTVoronoi3>(avoro, "voronoi3-dt.off");

	t.start();
	avoro.omp_compute(dt, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyDTVoronoi3>(avoro, "voronoi3-dt-omp.off");

	xlog("2 knn voronoi:");

	MyKNN3 knn;
	MyKnnVoronoi3 bvoro;

	t.start();
	knn.set_points(&sites[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	bvoro.set_minneigh(32);

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyKnnVoronoi3>(bvoro, "voronoi3-knn.off");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyKnnVoronoi3>(bvoro, "voronoi3-knn-omp.off");

	return 0;
}
