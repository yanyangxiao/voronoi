
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "trimesh.h"
#include "voronoi/lift.h"
#include "voronoi/facet.hpp"

// 1 rt method
#include "CGAL-rt3.h"
#include "voronoi/voronoi-dt.hpp"
// 2 knn method
#include "knn.hpp"
#include "voronoi/power-knn.hpp"

#include "timer.h"
#include "xlog.h"

int main(int argc, char *argv[])
{
	typedef double                 MyFloat;
	typedef xyy::Facet<3, MyFloat> MyFacet3;

	typedef MyCGALRT3                                      MyRT3;
	typedef xyy::DelaunayVoronoi<MyRT3, MyFacet3, MyFloat> MyRTPowers;

	typedef xyy::KNN<4, MyFloat>                           MyKNN4;
	typedef xyy::KnnPower<MyKNN4, MyFacet3, MyFloat>       MyKnnPowers;

	MyTrimesh tmesh;
	if (!OpenMesh::IO::read_mesh(tmesh, argv[1]))
		return 0;

	std::vector<MyFacet3> domain;
	convert_mesh<MyTrimesh, MyFacet3, MyFloat>(tmesh, domain);

	std::vector<MyFloat> sites;
	load_sites<4, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 4;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 rt power:");

	MyRT3 rt;
	MyRTPowers avoro;

	Timer t;
	t.start();
	rt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(rt, domain);
	t.stop();
	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyRTPowers>(avoro, "powers-rt.off");

	t.start();
	avoro.omp_compute(rt, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyRTPowers>(avoro, "powers-rt-omp.off");

	xlog("2 knn power:");

	MyKNN4 knn;
	MyKnnPowers bvoro;

	// lift
	std::vector<MyFloat> lifts;
	lift_sites<4, MyFloat>(sites, lifts);

	t.start();
	knn.set_points(&lifts[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyKnnPowers>(bvoro, "powers-knn.off");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_voronois<MyFloat, MyKnnPowers>(bvoro, "powers-knn-omp.off");

	return 0;
}
