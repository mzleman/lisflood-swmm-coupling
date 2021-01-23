#pragma once
#include <netcdf.h>
#include "lisflood.h"
#include "nc_tools.h"
#include "consts.h"



template<class Allocator = std::allocator<NUMERIC_TYPE>>
class DynamicRain
{
public:
	DynamicRain
	(
		const char* filename,
		int verbose,
		const Allocator& = Allocator()
	);

	/**
	 * Check rainfall data has same origin as the DEM.
	 */
	bool has_same_origin
	(
		Pars* Parptr
	);

	/**
	 * Check rainfall tile dimensions are an integer multiple of
	 * the DEM cell dimensions.
	 */
	bool is_tile_size_multiple_of_grid
	(
		Pars* Parptr
	);

	/**
	 * Update the current solver time, which loads new rainfall
	 * data if necessary.
	 */
	void update_time
	(
		NUMERIC_TYPE t
	);

	/**
		 * Update the current solver time, which loads new rainfall
		 * data if necessary. For GGM returns the bool yes if updated
		 */
		 //bool update_time_SGC
		 //(
		 //	NUMERIC_TYPE t
		 //);

		 /**
		  * The rainfall rate (m/s) at cell (i, j).
		  * If the cell is outside the rainfall domain, the rainfall rate is zero.
		  */
	NUMERIC_TYPE rate_at_cell
	(
		Pars* Parptr,
		int i,
		int j
	);

	NUMERIC_TYPE rate_at_cell
	(
		double cell_x,
		double cell_y
	);


	//NUMERIC_TYPE rate_at_cell_SGC
	//(
	//	int i,
	//	int j,
	//	const NUMERIC_TYPE dx_col,
	//	const NUMERIC_TYPE dy_col,
	//	NUMERIC_TYPE tly
	//);

	void update_H
	(
		Pars* Parptr,
		Solver *Solverptr,
		Arrays *Arrptr
	);

	//void update_rain_grid_SGM
	//(
	//	const NUMERIC_TYPE curr_time,
	//	NUMERIC_TYPE* rain_grid,
	//	const NUMERIC_TYPE* dem_grid,
	//	WetDryRowBound* wet_dry_bounds,
	//	const NUMERIC_TYPE* dx_col,
	//	const NUMERIC_TYPE* dy_col,
	//	const NUMERIC_TYPE tyl,
	//	const int grid_rows,
	//	const int grid_cols_padded,
	//	const NUMERIC_TYPE* cell_area_col
	//);

	inline bool enabled()
	{
		return enabled_;
	}

	inline const Geometry& geometry()
	{
		return geometry_;
	}

	/**
	 * Spatial map of rainfall rate (m/s) for the current time.
	 */
	inline NUMERIC_TYPE* data()
	{
		return netcdf_.data;
	}

	~DynamicRain();

private:
//public:
	bool enabled_;
	NetCDFVariable netcdf_;
	Geometry geometry_;
	Allocator allocator_;
};



template<class Allocator>
DynamicRain<Allocator>::DynamicRain
(
	const char* filename,
	int verbose,
	const Allocator& allocator
)
	:allocator_(allocator)
{
	enabled_ = (strlen(filename) > 0);
	if (!enabled_) return; read_file_netCDF_start(filename, "rainfall_depth", &netcdf_);
	netcdf_.data = allocator_.allocate(netcdf_.xlen * netcdf_.ylen);
	geometry_.xsz = netcdf_.xlen;
	geometry_.ysz = netcdf_.ylen;
	geometry_.dx = netcdf_.xs[1] - netcdf_.xs[0];
	geometry_.dy = netcdf_.ys[0] - netcdf_.ys[1];
	geometry_.blx = netcdf_.xs[0] - geometry_.dx / C(2.0);
	geometry_.tly = netcdf_.ys[0] + geometry_.dy / C(2.0);
	geometry_.bly = geometry_.tly - geometry_.ysz*geometry_.dy;

	if (verbose == ON) {
		printf("\nDynamic Rain Initialized.\n");
		cout << "blx:\t" << geometry_.blx << endl;
		cout << "bly:\t" << geometry_.bly << endl;
		cout << "dx:\t" << geometry_.dx << endl;
		cout << "dy:\t" << geometry_.dy << endl;
		cout << "X_size:" << geometry_.xsz << endl;
		cout << "Y_size:" << geometry_.ysz << endl;
		cout << "time series:" << endl;
		for (int i = 0; i < netcdf_.tlen; i++) {
			printf("%6.2f\t", netcdf_.times[i]);
		}
		printf("\n");
	}

	update_time(C(0.0));

}


template<class Allocator>
bool DynamicRain<Allocator>::has_same_origin
(
	Pars* Parptr
)
{
	return std::abs(Parptr->blx - geometry_.blx) < EPSILON &&
		std::abs(Parptr->bly - geometry_.bly) < EPSILON;
}

template<class Allocator>
bool DynamicRain<Allocator>::is_tile_size_multiple_of_grid
(
	Pars* Parptr
)
{
	NUMERIC_TYPE dx_ratio = geometry_.dx / Parptr->dx;
	NUMERIC_TYPE dy_ratio = geometry_.dy / Parptr->dy;

	return std::abs(std::round(dx_ratio) - dx_ratio) <= EPSILON &&
		std::abs(std::round(dy_ratio) - dy_ratio) <= EPSILON;
}

template<class Allocator>
void DynamicRain<Allocator>::update_time
(
	NUMERIC_TYPE t
)
{
	if (!enabled_) return;

	if (read_file_netCDF(&netcdf_, t / C(3600.0) /* s to hr */))
	{
		netcdf_.dt *= C(3600.0); // hr to s

		// convert from rainfall accumulated over dt (mm) to rainfall rate (m/s)
		for (int i = 0; i < geometry_.xsz*geometry_.ysz; i++)
		{
			netcdf_.data[i] /= netcdf_.dt; // £¿ ÀÛ»ý½µÓêÁ¿£¿
			netcdf_.data[i] /= C(1000.0); // mm to m
		}
	}
}


template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell
(
	Pars* Parptr,
	int i,
	int j
)
{
	if (!enabled_) return C(0.0);

	int tile_i = i * Parptr->dx / geometry_.dx;

	NUMERIC_TYPE top_gap = Parptr->tly - geometry_.tly;
	int tile_j = (j - top_gap / Parptr->dy) * Parptr->dy / geometry_.dy;

	if (tile_i < geometry_.xsz && tile_j >= 0)
	{
		return netcdf_.data[tile_j*geometry_.xsz + tile_i];
	}
	else
	{
		return C(0.0);
	}
}

template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell
(
	double cell_x,
	double cell_y
)
{
	if (!enabled_) return C(0.0);

	int tile_i = int( (cell_x - geometry_.blx) / geometry_.dx ),
		tile_j = int( (geometry_.tly - cell_y) / geometry_.dy );

	if (tile_i < 0 || tile_i >= netcdf_.xlen || tile_j < 0 || tile_j > netcdf_.ylen) {
		return C(0.0);
	}
	else {
		return netcdf_.data[tile_j*geometry_.xsz + tile_i];
	}

}


template<class Allocator>
void DynamicRain<Allocator>::update_H
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
#if _NETCDF == 1
	if (!enabled_) return;

	double x, y;
	int index;
	update_time(Solverptr->t);
	NUMERIC_TYPE total_rain_mass = C(0.0);
#pragma omp parallel for num_threads(Solverptr->ThreadNum) reduction (+:total_rain_mass) private(x, y, index) 
	for (int j = 0; j < Parptr->ysz; j++)
	{
		for (int i = 0; i < Parptr->xsz; i++)
		{
			index = j * Parptr->xsz + i;
			x = Arrptr->X_Coordinates[index];
			y = Arrptr->Y_Coordinates[index];
			if (fabs(x - NULLVAL) < EPSILON || fabs(y - NULLVAL) < EPSILON) continue;
			//NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
			//if (fabs(Z - Parptr->nodata_elevation) < 1e-9) continue;

			NUMERIC_TYPE& H = Arrptr->H[index];
			//NUMERIC_TYPE inc = rate_at_cell(Parptr, i, j) * Solverptr->Tstep;
			NUMERIC_TYPE inc = rate_at_cell(x, y) * Solverptr->Tstep;
			H += inc;
			total_rain_mass += inc * Parptr->dA;
		}
	}
#endif
	Parptr->RainTotalLoss += total_rain_mass;
}

template<class Allocator>
DynamicRain<Allocator>::~DynamicRain()
{
	if (!enabled_) return;
	free(netcdf_.xs);
	free(netcdf_.ys);
	free(netcdf_.times);
	allocator_.deallocate(netcdf_.data, netcdf_.xlen * netcdf_.ylen);
	CloseNetCDF(netcdf_.ncid);
}


