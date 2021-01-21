#pragma once
#include "lisflood.h"


#if _NUMERIC_MODE == 1
#define NC_GET_VARA_NUMERIC_TYPE nc_get_vara_double
#else
#define NC_GET_VARA_NUMERIC_TYPE nc_get_vara_float
#endif

void CloseNetCDF(int ncid);
bool read_file_netCDF(NetCDFVariable *ncvar, NUMERIC_TYPE curr_time);
void read_file_netCDF_start(const char *ncfilename, const char *varname, NetCDFVariable *ncvar);
