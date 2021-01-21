
#include "global.h"
#include "dynamic_rain.h"
#include "consts.h"

void initCoordinates(Pars *Parptr, Arrays *Arrptr, int verbose) {
	if (verbose == ON)  printf("initialize X,Y coordinates array.");
	double dx = Parptr->dx,
		dy = Parptr->dy,
		blx = Parptr->blx,
		tly = Parptr->tly;

	int size_x = Parptr->xsz,
		size_y = Parptr->ysz,
		count = 0, index;

	Arrptr->X_Coordinates = new double[size_x * size_y]();
	Arrptr->Y_Coordinates = new double[size_x * size_y]();

	for (int i = 0; i < size_y; i++) {
		for (int j = 0; j < size_x; j++) {
			index = j + i * size_x;
			if (fabs(Arrptr->DEM[index] - INFINIT) < EPSILON) {
				Arrptr->X_Coordinates[index] = NULLVAL;
				Arrptr->Y_Coordinates[index] = NULLVAL;
			}
			else {
				Arrptr->X_Coordinates[index] = blx + 0.5 * dx + j * dx;
				Arrptr->Y_Coordinates[index] = tly - 0.5 * dy - i * dy;
				++count;
			}
		}
	}
	if (verbose == ON)  printf("\nDone. %d cells.\n", count);
}