
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "voronoi/lift.h"
#include "voronoi/facet.hpp"

// rt method
#include "CGAL-rt2.h"
#include "voronoi/voronoi-dt.hpp"
// knn method
#include "knn.hpp"
#include "voronoi/power-knn.hpp"

#include "timer.h"
#include "xlog.h"

template <int SiteDIM, typename REAL, typename Voro2>
bool save_power2(const REAL *sites, int vnb, const Voro2 &voro, const char *svgname)
{
	FILE *svg = fopen(svgname, "wb");
	if (!svg)
		return false;

	REAL bbox[4] = { 1e10, -1e10, 1e10, -1e10 };
	for (int v = 0; v < (int)voro.cells_number(); ++v)
	{
		const Voro2::Cell& cell = voro.cell(v);
		for (int f = 0; f < (int)cell.size(); ++f)
		{
			for (int k = 0; k < (int)cell[f].points_number(); ++k)
			{
				const REAL *p = cell[f].point(k);
				if (p[0] < bbox[0])
					bbox[0] = p[0];
				if (p[0] > bbox[1])
					bbox[1] = p[0];
				if (p[1] < bbox[2])
					bbox[2] = p[1];
				if (p[1] > bbox[3])
					bbox[3] = p[1];
			}
		}
	}

	REAL dx = bbox[1] - bbox[0];
	REAL dy = bbox[3] - bbox[2];

	REAL scale = 1000.0 / dx;

	int margin = 100;
	int width = int(dx * scale) + margin;
	int height = int(dy * scale) + margin;

	fprintf(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n  width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
		width, height, width, height);
	fprintf(svg, "<metadata>\n");
	fprintf(svg, "  Created by Xiao, Yanyang (yanyangxiaoxyy@gmail.com)\n");
	fprintf(svg, "  All copyrights reserved\n");
	fprintf(svg, "</metadata>\n");

	fprintf(svg, "<g id=\"cells\">\n");
	for (int v = 0; v < (int)voro.cells_number(); ++v)
	{
		const Voro2::Cell& cell = voro.cell(v);

		int r = rand() % 200 + 50;
		int g = rand() % 200 + 50;
		int b = rand() % 200 + 50;

		fprintf(svg, "<g fill=\"rgb(%d, %d, %d)\">\n", r, g, b);

		for (int f = 0; f < (int)cell.size(); ++f)
		{
			fprintf(svg, "<path d=\"M ");
			for (int k = 0; k < (int)cell[f].points_number(); ++k)
			{
				const REAL *p = cell[f].point(k);

				REAL svgx = (p[0] - bbox[0]) * scale + margin / 2;
				REAL svgy = height - ((p[1] - bbox[2]) * scale + margin / 2);

				fprintf(svg, "%f %f ", svgx, svgy);
				if (k < (int)cell[f].points_number() - 1)
					fprintf(svg, "L ");
			}
			fprintf(svg, " Z\"/>\n");
		}

		fprintf(svg, "</g>\n");
	}
	fprintf(svg, "</g>\n");

	fprintf(svg, "<g id=\"edges\" stroke=\"rgb(0, 0, 0)\" stroke-width=\"3\">\n");
	for (int v = 0; v < (int)voro.cells_number(); ++v)
	{
		const Voro2::Cell& cell = voro.cell(v);

		for (int f = 0; f < (int)cell.size(); ++f)
		{
			for (int k = 0; k < (int)cell[f].points_number(); ++k)
			{
				const REAL *p = cell[f].point(k);
				int sn = cell[f].wall_siteneigh(k);
				int dn = cell[f].wall_domneigh(k);
				if ((dn < -1) || (sn < v && sn > -1))
					continue;

				int next = (k + 1) % (int)cell[f].points_number();
				const REAL *q = cell[f].point(next);

				REAL svgx1 = (p[0] - bbox[0]) * scale + margin / 2;
				REAL svgy1 = height - ((p[1] - bbox[2]) * scale + margin / 2);
				REAL svgx2 = (q[0] - bbox[0]) * scale + margin / 2;
				REAL svgy2 = height - ((q[1] - bbox[2]) * scale + margin / 2);

				fprintf(svg, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\"/>", svgx1, svgy1, svgx2, svgy2);
			}
		}
	}
	fprintf(svg, "</g>\n");

	fprintf(svg, "<g id=\"negative sites\" fill=\"rgb(50, 50, 250)\">\n");
	for (int v = 0; v < vnb; ++v)
	{
		const REAL *vpos = &sites[SiteDIM * v];
		if (vpos[SiteDIM - 1] >= 0)
			continue;

		REAL svgx = (vpos[0] - bbox[0]) * scale + margin / 2;
		REAL svgy = height - ((vpos[1] - bbox[2]) * scale + margin / 2);
		fprintf(svg, "<circle cx=\"%f\" cy=\"%f\" r=\"4\"/>\n", svgx, svgy);
	}
	fprintf(svg, "</g>\n");

	fprintf(svg, "<g id=\"positive sites\" fill=\"rgb(250, 50, 50)\">\n");
	for (int v = 0; v < vnb; ++v)
	{
		const REAL *vpos = &sites[SiteDIM * v];
		if (vpos[SiteDIM - 1] < 0)
			continue;

		REAL svgx = (vpos[0] - bbox[0]) * scale + margin / 2;
		REAL svgy = height - ((vpos[1] - bbox[2]) * scale + margin / 2);
		fprintf(svg, "<circle cx=\"%f\" cy=\"%f\" r=\"4\"/>\n", svgx, svgy);
	}
	fprintf(svg, "</g>\n");

	fprintf(svg, "</svg>\n");
	fclose(svg);

	return true;
}

int main(int argc, char *argv[])
{
	typedef double                    MyFloat;
	typedef xyy::Facet<2, MyFloat>    MyFacet2;

	typedef MyCGALRT2                                      MyRT2;
	typedef xyy::DelaunayVoronoi<MyRT2, MyFacet2, MyFloat> MyRTPower2;

	typedef xyy::KNN<3, MyFloat>                           MyKNN3;
	typedef xyy::KnnPower<MyKNN3, MyFacet2, MyFloat>       MyKnnPower2;

	std::vector<MyFacet2> domain;
	load_domain<MyFacet2, MyFloat>(argv[1], domain);
	
	std::vector<MyFloat> sites;
	load_sites<3, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 3;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 rt power:");

	MyRT2 rt;
	MyRTPower2 avoro;

	Timer t;
	t.start();
	rt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(rt, domain);
	t.stop();
	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_power2<3, MyFloat, MyRTPower2>(&sites[0], vnb, avoro, "power2-rt.svg");

	t.start();
	avoro.omp_compute(rt, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_power2<3, MyFloat, MyRTPower2>(&sites[0], vnb, avoro, "power2-rt-omp.svg");

	xlog("2 knn power:");

	MyKNN3 knn;
	MyKnnPower2 bvoro;

	// lift
	std::vector<MyFloat> lifts;
	lift_sites<3, MyFloat>(sites, lifts);

	t.start();
	knn.set_points(&lifts[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute power: %.10f seconds", t.elapsed_time());

	save_power2<3, MyFloat, MyKnnPower2>(&sites[0], vnb, bvoro, "power2-knn.svg");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute power: %.10f seconds", t.elapsed_time());

	save_power2<3, MyFloat, MyKnnPower2>(&sites[0], vnb, bvoro, "power2-knn-omp.svg");

	return 0;
}
