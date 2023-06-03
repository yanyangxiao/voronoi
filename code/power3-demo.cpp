
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "voronoi/lift.h"
#include "voronoi/polymesh.hpp"

// dt
#include "CGAL-rt3.h"
#include "voronoi/voronoi-dt.hpp"
// knn
#include "knn.hpp"
#include "voronoi/power-knn.hpp"

#include "timer.h"
#include "xlog.h"

int main(int argc, char *argv[])
{
	typedef double                    MyFloat;
	typedef xyy::Polymesh<MyFloat>    MyPolymesh;

	typedef MyCGALRT3                                        MyRT3;
	typedef xyy::DelaunayVoronoi<MyRT3, MyPolymesh, MyFloat> MyRTPower3;

	typedef xyy::KNN<4, MyFloat>                             MyKNN4;
	typedef xyy::KnnPower<MyKNN4, MyPolymesh, MyFloat>       MyKnnPower3;

	std::vector<MyPolymesh> domain;
	load_tets<MyPolymesh, MyFloat>(argv[1], domain);

	std::vector<MyFloat> sites;
	load_sites<4, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 4;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 rt power:");

	MyRT3 dt;
	MyRTPower3 avoro;

	Timer t;
	t.start();
	dt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(dt, domain);
	t.stop();
	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyRTPower3>(avoro, "power3-rt.off");

	t.start();
	avoro.omp_compute(dt, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyRTPower3>(avoro, "power3-rt-omp.off");

	xlog("2 knn power:");

	MyKNN4 knn;
	MyKnnPower3 bvoro;

	// lift
	std::vector<MyFloat> lifts;
	lift_sites<4, MyFloat>(sites, lifts);

	t.start();
	knn.set_points(&lifts[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	bvoro.set_minneigh(32);

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyKnnPower3>(bvoro, "power3-knn.off");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_voronoi3<MyFloat, MyKnnPower3>(bvoro, "power3-knn-omp.off");

	return 0;
}
