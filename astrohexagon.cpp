/************************************************************************
astrohexagon.cpp
Set up a hexagonal mesh of astrocytes in a specified plane
Provides input to cone detector (ICeD.cpp)
TWS May 2021
**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

//sets up 7-hexagon mesh of astrocytes centered on origin
void astrosetup(void) {
	extern float INLzonetop, INLzonebottom, ahex;		
	extern float ** astropts;
	int i;
	float cc = sqrt(3.) / 2.;
	float ahexp = ahex * cc, ahex2 = ahex / 2.;

	astropts[1][1] = ahexp;	//start of segment
	astropts[2][1] = ahex2;
	astropts[4][1] = 0.;	//end of segment
	astropts[5][1] = ahex;
	astropts[1][7] = ahexp;	//start of segment
	astropts[2][7] = ahex2;
	astropts[4][7] = 2. * ahexp;	//end of segment
	astropts[5][7] = ahex;
	astropts[1][13] = 2. * ahexp;	//start of segment
	astropts[2][13] = ahex;
	astropts[4][13] = 2. * ahexp;	//end of segment
	astropts[5][13] = 2. * ahex;
	astropts[1][19] = 2. * ahexp;	//start of segment
	astropts[2][19] = 2. * ahex;
	astropts[4][19] = ahexp;	//end of segment
	astropts[5][19] = 5. * ahex2;
	astropts[1][25] = ahexp;	//start of segment
	astropts[2][25] = 5. * ahex2;
	astropts[4][25] = 0.;	//end of segment
	astropts[5][25] = 2. * ahex;
	for (i = 1; i <= 30; i++) {
		if (i % 6 != 1) {			//rotate through pi/3
			astropts[1][i] = astropts[1][i - 1] * 0.5 - astropts[2][i - 1] * cc;
			astropts[2][i] = astropts[1][i - 1] * cc + astropts[2][i - 1] * 0.5;
			astropts[4][i] = astropts[4][i - 1] * 0.5 - astropts[5][i - 1] * cc;
			astropts[5][i] = astropts[4][i - 1] * cc + astropts[5][i - 1] * 0.5;
		}
		astropts[3][i] = INLzonetop;
		astropts[6][i] = INLzonetop;
		astropts[1][i + 30] = astropts[1][i];
		astropts[2][i + 30] = astropts[2][i];
		astropts[3][i + 30] = INLzonebottom;
		astropts[4][i + 30] = astropts[4][i];
		astropts[5][i + 30] = astropts[5][i];
		astropts[6][i + 30] = INLzonebottom;
	}
	//FILE* ofp;		//for testing
	//ofp = fopen("test.ps", "w");
	//fprintf(ofp, "%%!PS-Adobe-2.0\n");
	//fprintf(ofp, "/m {moveto} def\n");
	//fprintf(ofp, "/l {lineto} def\n");
	//fprintf(ofp, "/s {stroke} def\n");
	//float xmin = 0., ymin = 0., scalefac = 1., cx = 200., cy = 200.;
	//fprintf(ofp, "/mx {%g sub %g mul %g add} def\n", xmin, scalefac, cx);
	//fprintf(ofp, "/my {%g sub %g mul %g add} def\n", ymin, scalefac, cy);
	//for (i = 1; i <= 30; i++)
	//	fprintf(ofp, "%g mx %g my m %g mx %g my l s\n", astropts[1][i], astropts[2][i], astropts[4][i], astropts[5][i]);
	//fclose(ofp);
}

void astrohexagon(float x, float y, float *xc, float *yc)
{
	extern float ahex, alx, aly;
	float ahexp = ahex * sqrt(3.) / 2.;
	float x1, y1, fx, fy, ixc, iyc, x1f, y1f;

	//locate center of hexagon containing (x,y) 
	//xcenter = ixc * ahexp, ycenter = iyc * ahex * 3/2
	x1 = (x - alx / 2.)/ 2. / ahexp;
	x1f = floor(x1);
	fx = x1 - x1f;
	y1 = (y - aly / 2.)/ 3. / ahex;
	y1f = floor(y1);
	fy = y1 - y1f;
	if (fx <= 0.5) {
		if (fx + 3. * fy <= 1.) {
			ixc = 2 * x1f;
			iyc = 2 * y1f;
		}
		else if (-fx + 3. * fy <= 2.) {
			ixc = 2 * x1f + 1;
			iyc = 2 * y1f + 1;
		}
		else {
			ixc = 2 * x1f;
			iyc = 2 * y1f + 2;
		}
	}
	else {
		if (-fx + 3. * fy <= 0.) {
			ixc = 2 * x1f + 2;
			iyc = 2 * y1f;
		}
		else if (fx + 3. * fy <= 3.) {
			ixc = 2 * x1f + 1;
			iyc = 2 * y1f + 1;
		}
		else {
			ixc = 2 * x1f + 2;
			iyc = 2 * y1f + 2;
		}
	}

	*xc = ixc * ahexp + alx / 2.;
	*yc = iyc * ahex * 1.5 + aly / 2.;

	float dist = sqrt(SQR(*xc - x) + SQR(*yc - y));
	if (dist > ahex) printf("***Error in astrohexagon\n");
}
