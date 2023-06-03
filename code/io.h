
#ifndef IO_H
#define IO_H

#include <vector>

template <typename Facet, typename REAL>
bool load_domain(const char *filename, std::vector<Facet> &domain)
{
	domain.clear();

	std::ifstream file(filename);
	if (!file)
		return false;

	int vnb, fnb;
	file >> vnb >> fnb;

	int dim = Facet::Dimension;

	std::vector<REAL> points;
	points.resize(dim * vnb);

	for (int k = 0; k < vnb; ++k)
	{
		for (int d = 0; d < dim; ++d)
			file >> points[dim * k + d];
	}

	domain.resize(fnb);

	for (int f = 0; f < fnb; ++f)
	{
		int fvnb;
		file >> fvnb;

		for (int i = 0; i < fvnb; ++i)
		{
			int fv, fneigh;
			file >> fv >> fneigh;

			domain[f].add_point(&points[dim * fv], fneigh, -1);
		}
	}

	file.close();

	return true;
}

template <typename TMesh, typename Facet, typename REAL>
void convert_mesh(const TMesh &tmesh, std::vector<Facet> &domain)
{
	domain.clear();
	domain.resize(tmesh.n_faces());

	for (int f = 0; f < (int)tmesh.n_faces(); ++f)
	{
		TMesh::FaceHandle fh = tmesh.face_handle(f);
		TMesh::HalfedgeHandle fhe = tmesh.halfedge_handle(fh);
		TMesh::VertexHandle vh = tmesh.from_vertex_handle(fhe);
		TMesh::Point vp = tmesh.point(vh);
		int fneigh = tmesh.opposite_face_handle(fhe).idx();

		domain[f].add_point(&vp[0], fneigh, -1);

		fhe = tmesh.next_halfedge_handle(fhe);
		vh = tmesh.from_vertex_handle(fhe);
		vp = tmesh.point(vh);
		fneigh = tmesh.opposite_face_handle(fhe).idx();

		domain[f].add_point(&vp[0], fneigh, -1);

		fhe = tmesh.next_halfedge_handle(fhe);
		vh = tmesh.from_vertex_handle(fhe);
		vp = tmesh.point(vh);
		fneigh = tmesh.opposite_face_handle(fhe).idx();

		domain[f].add_point(&vp[0], fneigh, -1);
	}
}

template <typename SingleCell, typename REAL>
bool load_tets(const char *filename, std::vector<SingleCell> &domain)
{
	domain.clear();

	std::ifstream file(filename);
	if (!file)
		return false;

	int pnb, cnb;
	file >> pnb >> cnb;

	std::vector<REAL> points(3 * pnb);
	for (int v = 0; v < pnb; ++v)
	{
		file >> points[3 * v] >> points[3 * v + 1] >> points[3 * v + 2];
	}

	std::vector<std::vector<int>> tets;
	std::vector<std::vector<int>> neighs;

	tets.resize(cnb);
	neighs.resize(cnb);

	for (int c = 0; c < cnb; ++c)
	{
		tets[c].clear();
		tets[c].resize(4);

		file >> tets[c][0] >> tets[c][1] >> tets[c][2] >> tets[c][3];
	}

	for (int c = 0; c < cnb; ++c)
	{
		neighs[c].clear();
		neighs[c].resize(4);

		file >> neighs[c][0] >> neighs[c][1] >> neighs[c][2] >> neighs[c][3];
	}

	file.close();

	domain.resize(cnb);

	for (int c = 0; c < cnb; ++c)
	{
		int va = tets[c][0];
		int vb = tets[c][1];
		int vc = tets[c][2];
		int vd = tets[c][3];

		int na = neighs[c][0];
		int nb = neighs[c][1];
		int nc = neighs[c][2];
		int nd = neighs[c][3];

		domain[c].begin_facet();
		domain[c].add_point(&points[3 * va]);
		domain[c].add_point(&points[3 * vc]);
		domain[c].add_point(&points[3 * vb]);
		domain[c].end_facet(nd, -1);

		domain[c].begin_facet();
		domain[c].add_point(&points[3 * va]);
		domain[c].add_point(&points[3 * vb]);
		domain[c].add_point(&points[3 * vd]);
		domain[c].end_facet(nc, -1);

		domain[c].begin_facet();
		domain[c].add_point(&points[3 * va]);
		domain[c].add_point(&points[3 * vd]);
		domain[c].add_point(&points[3 * vc]);
		domain[c].end_facet(nb, -1);

		domain[c].begin_facet();
		domain[c].add_point(&points[3 * vb]);
		domain[c].add_point(&points[3 * vc]);
		domain[c].add_point(&points[3 * vd]);
		domain[c].end_facet(na, -1);
	}

	return true;
}

template <int DIM, typename REAL>
bool load_sites(const char *filename, std::vector<REAL> &sites)
{
	std::ifstream file(filename);
	if (!file)
		return false;

	int vnb;
	file >> vnb;

	sites.clear();
	sites.resize(DIM * vnb);
	for (int v = 0; v < vnb; ++v)
	{
		for (int d = 0; d < DIM; ++d)
			file >> sites[DIM * v + d];
	}

	file.close();

	return true;
}

template <typename REAL, typename Voros>
bool save_voronois(const Voros &voro, const char *offname)
{
	std::ofstream file(offname);
	if (!file)
		return false;

	int pnb = 0;
	int fnb = 0;

	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voros::Cell &cell = voro.cell(v);

		for (int f = 0; f < (int)cell.size(); ++f)
			pnb += (int)cell[f].points_number();

		fnb += (int)cell.size();
	}

	file << "OFF" << std::endl;
	file << pnb << " " << fnb << " " << 0 << std::endl;

	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voros::Cell &cell = voro.cell(v);
		for (int f = 0; f < (int)cell.size(); ++f)
		{
			for (int i = 0; i < (int)cell[f].points_number(); ++i)
			{
				const REAL *p = cell[f].point(i);
				file << p[0] << " " << p[1] << " " << p[2] << std::endl;
			}
		}
	}

	int index = 0;
	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voros::Cell &cell = voro.cell(v);

		int r = rand() % 180 + 50;
		int g = rand() % 180 + 50;
		int b = rand() % 180 + 50;

		for (int f = 0; f < (int)cell.size(); ++f)
		{
			file << cell[f].points_number();

			for (int i = 0; i < (int)cell[f].points_number(); ++i)
			{
				file << " " << i + index;
			}

			file << " " << float(r) / 255.0f << " " << float(g) / 255.0 << " " << float(b) / 255.0f << std::endl;

			index += (int)cell[f].points_number();
		}
	}

	file.close();

	return true;
}

template <typename REAL, typename Voro3>
bool save_voronoi3(const Voro3 &voro, const char *offname, double shrink = 0.9)
{
	std::ofstream file(offname);
	if (!file)
		return false;

	int pnb = 0;
	int fnb = 0;

	std::vector<REAL> cellcents;
	cellcents.resize(3 * voro.cells_number());

	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voro3::Cell &cell = voro.cell(v);
		if (cell.empty())
			continue;

		cellcents[3 * v] = 0;
		cellcents[3 * v + 1] = 0;
		cellcents[3 * v + 2] = 0;

		for (int it = 0; it < (int)cell.size(); ++it)
		{
			for (int f = 0; f < cell[it].walls_number(); ++f)
			{
				int sn = cell[it].wall_siteneigh(f);
				int dn = cell[it].wall_domneigh(f);

				if (sn < 0 && dn > -1)
					continue;

				pnb += (int)cell[it].facet_size(f);
				++fnb;
			}

			REAL cent[3];
			cell[it].center(cent);

			cellcents[3 * v] += cent[0];
			cellcents[3 * v + 1] += cent[1];
			cellcents[3 * v + 2] += cent[2];
		}

		cellcents[3 * v] /= (int)cell.size();
		cellcents[3 * v + 1] /= (int)cell.size();
		cellcents[3 * v + 2] /= (int)cell.size();
	}

	file << "OFF" << std::endl;
	file << pnb << " " << fnb << " " << 0 << std::endl;

	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voro3::Cell &cell = voro.cell(v);

		for (int it = 0; it < (int)cell.size(); ++it)
		{
			for (int f = 0; f < cell[it].walls_number(); ++f)
			{
				int sn = cell[it].wall_siteneigh(f);
				int dn = cell[it].wall_domneigh(f);

				if (sn < 0 && dn > -1)
					continue;

				for (int i = cell[it].facet_begin(f); i < cell[it].facet_end(f); ++i)
				{
					const REAL *p = cell[it].point(i);
					file << p[0] * shrink + cellcents[3 * v] * (1.0 - shrink) << " "
						<< p[1] * shrink + cellcents[3 * v + 1] * (1.0 - shrink) << " "
						<< p[2] * shrink + cellcents[3 * v + 2] * (1.0 - shrink) << std::endl;
				}
			}
		}
	}

	int index = 0;
	for (int v = 0; v < voro.cells_number(); ++v)
	{
		const Voro3::Cell &cell = voro.cell(v);

		int r = rand() % 180 + 50;
		int g = rand() % 180 + 50;
		int b = rand() % 180 + 50;

		for (int it = 0; it < (int)cell.size(); ++it)
		{
			for (int f = 0; f < cell[it].walls_number(); ++f)
			{
				int sn = cell[it].wall_siteneigh(f);
				int dn = cell[it].wall_domneigh(f);

				if (sn < 0 && dn > -1)
					continue;

				file << cell[it].facet_size(f);

				for (int fv = cell[it].facet_begin(f); fv < cell[it].facet_end(f); ++fv)
				{
					file << " " << fv - cell[it].facet_begin(f) + index;
				}

				file << " " << float(r) / 255.0f << " " << float(g) / 255.0 << " " << float(b) / 255.0f << std::endl;
			
				index += cell[it].facet_size(f);
			}
		}
	}

	file.close();

	return true;
}

#endif
