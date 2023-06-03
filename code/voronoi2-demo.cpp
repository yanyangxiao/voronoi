
#include <stdlib.h>
#include <fstream>

#include "io.h"
#include "voronoi/facet.hpp"

// dt method
#include "CGAL-dt2.h"
#include "voronoi/voronoi-dt.hpp"
// knn method
#include "knn.hpp"
#include "voronoi/voronoi-knn.hpp"

#include "timer.h"
#include "xlog.h"

template <typename REAL, typename Voro2>
bool save_voronoi2(const REAL *sites, int vnb, const Voro2 &voro, const char *svgname)
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

	fprintf(svg, "<g id=\"sites\" fill=\"rgb(250, 50, 50)\">\n");
	for (int v = 0; v < vnb; ++v)
	{
		const REAL *vpos = &sites[2 * v];
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

	typedef MyCGALDT2                                      MyDT2;
	typedef xyy::DelaunayVoronoi<MyDT2, MyFacet2, MyFloat> MyDTVoronoi2;

	typedef xyy::KNN<2, MyFloat>                           MyKNN2;
	typedef xyy::KnnVoronoi<MyKNN2, MyFacet2, MyFloat>     MyKnnVoronoi2;

	std::vector<MyFacet2> domain;
	load_domain<MyFacet2, MyFloat>(argv[1], domain);
	
	std::vector<MyFloat> sites;
	load_sites<2, MyFloat>(argv[2], sites);

	int vnb = (int)sites.size() / 2;

	xlog("input: regions = %d, sites = %d", (int)domain.size(), vnb);

	xlog("1 dt voronoi:");

	MyDT2 dt;
	MyDTVoronoi2 avoro;

	Timer t;
	t.start();
	dt.set_vertices(&sites[0], vnb);
	t.stop();
	xlog("compute triangulation: %.10f seconds", t.elapsed_time());

	t.start();
	avoro.compute(dt, domain);
	t.stop();
	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi2<MyFloat, MyDTVoronoi2>(&sites[0], vnb, avoro, "voronoi2-dt.svg");

	t.start();
	avoro.omp_compute(dt, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi2<MyFloat, MyDTVoronoi2>(&sites[0], vnb, avoro, "voronoi2-dt-omp.svg");

	xlog("2 knn voronoi:");

	MyKNN2 knn;
	MyKnnVoronoi2 bvoro;

	t.start();
	knn.set_points(&sites[0], vnb);
	t.stop();
	xlog("compute kdtree: %.10f seconds", t.elapsed_time());

	bvoro.set_minneigh(16);

	t.start();
	bvoro.compute(knn, domain);
	t.stop();

	xlog("compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi2<MyFloat, MyKnnVoronoi2>(&sites[0], vnb, bvoro, "voronoi2-knn.svg");

	t.start();
	bvoro.omp_compute(knn, domain);
	t.stop();

	xlog("omp_compute voronoi: %.10f seconds", t.elapsed_time());

	save_voronoi2<MyFloat, MyKnnVoronoi2>(&sites[0], vnb, bvoro, "voronoi2-knn-omp.svg");

	return 0;
}
