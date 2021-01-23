#include "nc_tools.h"
#include "global.h"

#if _NETCDF == 1
#include <netcdf.h>
#endif

/* Handle errors by printing an error message and exiting with a
* non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#if _NUMERIC_MODE == 1
#define nc_put_att_NUMERIC_TYPE nc_put_att_double
#define nc_put_var1_NUMERIC_TYPE nc_put_var1_double
#define NC_NUMERIC_TYPE NC_DOUBLE

#else
#define nc_put_att_NUMERIC_TYPE nc_put_att_float
#define nc_put_var1_NUMERIC_TYPE nc_put_var1_float
#define NC_NUMERIC_TYPE NC_FLOAT
#endif

#define UNITS_LABEL "units"
#define LONG_NAME_LABEL "long_name"
#define MISSING_VALUE_LABEL	"missing_value"

#define HEIGHT_UNITS_VALUE "meters" 
#define FLOW_UNITS_VALUE "meters3/second" 
#define VELOCITY_UNITS_VALUE "meters/second" 
#define TIME_DURATION_HOURS_UNITS_VALUE "hours"

#if _NETCDF == 1

void CloseNetCDF(int ncid)
{
	int retval;
	if (retval = nc_close(ncid)) ERR(retval);
}


void read_file_netCDF_start(const char *ncfilename, const char *varname,
	NetCDFVariable *ncvar)
{
	int recid, xid, yid;
	int retval;

	if ((retval = nc_open(ncfilename, NC_NOWRITE, &ncvar->ncid)))
		ERR(retval);
	if ((retval = nc_inq_varid(ncvar->ncid, varname, &ncvar->varid)))
		ERR(retval);
	if ((retval = nc_inq_varid(ncvar->ncid, "time", &recid)))
		ERR(retval);
	if ((retval = nc_inq_dimid(ncvar->ncid, "x", &xid)))
		ERR(retval);
	if ((retval = nc_inq_dimid(ncvar->ncid, "y", &yid)))
		ERR(retval);
	if ((retval = nc_inq_dimlen(ncvar->ncid, recid, &ncvar->tlen)))
		ERR(retval);
	if ((retval = nc_inq_dimlen(ncvar->ncid, xid, &ncvar->xlen)))
		ERR(retval);
	if ((retval = nc_inq_dimlen(ncvar->ncid, yid, &ncvar->ylen)))
		ERR(retval);

	size_t start[] = { 0 };
	size_t tcount[] = { ncvar->tlen };
	ncvar->times = (NUMERIC_TYPE*)calloc(ncvar->tlen, sizeof(NUMERIC_TYPE));
	if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, recid, start, tcount,
		ncvar->times))
		ERR(retval);
	ncvar->time_idx = -1;

	size_t xcount[] = { ncvar->xlen };
	ncvar->xs = (NUMERIC_TYPE*)calloc(ncvar->xlen, sizeof(NUMERIC_TYPE));
	if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, xid, start, xcount,
		ncvar->xs))
		ERR(retval);

	size_t ycount[] = { ncvar->ylen };
	ncvar->ys = (NUMERIC_TYPE*)calloc(ncvar->ylen, sizeof(NUMERIC_TYPE));
	if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, yid, start, ycount,
		ncvar->ys))
		ERR(retval);

}

// assume data are ordered left-to-right, top-to-bottom
// return true if new data needed to be read for the curr_time
bool read_file_netCDF(NetCDFVariable *ncvar, NUMERIC_TYPE curr_time)
{
	int retval;

	size_t old_idx = ncvar->time_idx;
	while (static_cast<int>(ncvar->time_idx) < static_cast<int>(ncvar->tlen) - 1 &&
		ncvar->times[ncvar->time_idx + 1] <= curr_time)
	{
		ncvar->time_idx++;
	}

	if (old_idx != ncvar->time_idx) {
		if (ncvar->time_idx < ncvar->tlen - 1)
		{
			ncvar->dt = ncvar->times[ncvar->time_idx + 1]
				- ncvar->times[ncvar->time_idx];
		}
		else  //
		{
			ncvar->dt = ncvar->times[ncvar->time_idx]
				- ncvar->times[ncvar->time_idx - 1];  // nc格式中降雨步长可以是不均匀的。
		}

		size_t start[] = { ncvar->time_idx, 0, 0 };
		size_t count[] = { 1, ncvar->ylen, ncvar->xlen };
		if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, ncvar->varid, start,
			count, ncvar->data)) ERR(retval);

		if (*verbose == ON) {
			printf("\n");
			cout << "Read netCDF rain file :" << endl;
			for (int i = 0; i < ncvar->ylen; i++) {
				for (int j = 0; j < ncvar->xlen; j++) {
					printf("%6.2f\t", ncvar->data[j + i * ncvar->xlen]);
				}
				printf("\n");
			}
			printf("\n");
		}		
	}

	return old_idx != ncvar->time_idx;
}

/*
netcdf_varid: must be defined in netcdf_state
*/
//void write_file_netCDF(NetCDFState * netcdf_state, int netcdf_varid,
//	int time_increment,
//	const NUMERIC_TYPE *data,
//	const int grid_cols, const int grid_rows, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size)
//{
//	int retval;
//	size_t count[3];
//	count[0] = 1; //time
//	count[1] = grid_rows; // y
//	count[2] = grid_cols; // x
//
//	size_t start[3];
//	start[0] = time_increment; //time
//	start[1] = 0; // y
//	start[2] = 0; // x
//
//	if ((retval = nc_put_vara(netcdf_state->ncid, netcdf_varid, start, count, data)))
//		ERR(retval);
//}

#endif
