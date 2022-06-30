/************************************************************************
analyzeresults for AngioAdapt
TWS, June 2020
With analysis of segment z-distribution
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void histogram(float* var, float* weight, int n, const char filename[]);

void analyzeresults()
{
	extern int imain,nnod, nseg, nnt, nnv, phaseseparation, * ista, * iend, * mainseg;
	extern int* nodtyp, * segtyp, ** nodseg, ** tisspoints, * noxy_z, mzz;
	extern float pi1, alz, INLzonetop, INLzonebottom, percentINLtop, percentINLbottom, percentneither;
	extern float* q, * qq, * diam, * lseg, * tau, * segpress, * hd, * dtmin, * sourcefac, * ds, * oxy_z;
	extern float** pt, ** pvseg, ** cnode;
	extern double* nodpress, * cond;

	int i, iseg, inod, maxdim, iz;
	float zval;
	float INLzonetopfrac, INLzonebottomfrac, NoINLzonefrac, totallength_weightfrac, totallength_weight;
	float* histogramdisplay, * histogramweight, * length_weight;
	char fname[80];
	FILE* foxy_z;

	maxdim = IMAX(nnod, nseg);
	maxdim = IMAX(maxdim, nnv);
	maxdim = IMAX(maxdim, nnt);
	histogramdisplay = vector(1, maxdim);
	histogramweight = vector(1, maxdim);
	length_weight = vector(1, nnod);

	//node weighting factors: 1/2 the sum of the lengths of the segments connected to the node
	for (inod = 1; inod <= nnod; inod++) {
		length_weight[inod] = 0;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) length_weight[inod] += 0.5 * lseg[iseg];
		}
	}

	//Node pressures
	i = 0;
	for (inod = 1; inod <= nnod; inod++) if (length_weight[inod] > 0.) {
		i++;
		histogramdisplay[i] = nodpress[inod];
		histogramweight[i] = length_weight[inod];
	}
	sprintf(fname, "Histogram-pressures%03i.out", imain);
	histogram(histogramdisplay, histogramweight, i, fname);

	//Segment flows
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(qq[iseg]);
		histogramweight[i] = lseg[iseg];
	}
	sprintf(fname, "Histogram-logflows%03i.out", imain);
	histogram(histogramdisplay, histogramweight, i, fname);

	//Segment velocities in mm/s
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(4000. / 60. * qq[iseg] / pi1 / SQR(diam[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	sprintf(fname, "Histogram-logvelocities%03i.out", imain);
	histogram(histogramdisplay, histogramweight, i, fname);

	//Shear stress
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && fabs(tau[iseg]) > 0.) {
		i++;
		histogramdisplay[i] = log10(fabs(tau[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	sprintf(fname, "Histogram-logstress%03i.out", imain);
	histogram(histogramdisplay, histogramweight, i, fname);

	//vessel PO2
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		i++;
		histogramdisplay[i] = pvseg[iseg][1];
		histogramweight[i] = lseg[iseg];
	}
	sprintf(fname, "Histogram-PO2v%03i.out", imain);
	histogram(histogramdisplay, histogramweight, i, fname);

	//tissue PO2
	for (i = 1; i <= nnt; i++) {
		histogramdisplay[i] = FMIN(pt[i][1], 99.9);
		histogramweight[i] = sourcefac[i];
	}
	sprintf(fname, "Histogram-PO2t%03i.out", imain);
	histogram(histogramdisplay, histogramweight, nnt, fname);

	//distance to nearest vessel
	for (i = 1; i <= nnt; i++) {
		histogramdisplay[i] = dtmin[i];
		histogramweight[i] = sourcefac[i];
	}
	sprintf(fname, "Histogram-Distance%03i.out", imain);
	histogram(histogramdisplay, histogramweight, nnt, fname);
	
	//hematocrit
	if (phaseseparation) {
		i = 0;
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			i++;
			histogramdisplay[i] = hd[iseg];
			histogramweight[i] = lseg[iseg];
		}
		sprintf(fname, "Histogram-hematocrit%03i.out", imain);
		histogram(histogramdisplay, histogramweight, i, fname);
	}

	// Z values
	INLzonetopfrac = 0.;
	INLzonebottomfrac = 0.;
	NoINLzonefrac = 0.;
	totallength_weightfrac = 0.;
	totallength_weight = 0;

	for (i = 1; i <= nnod; i++) {
		histogramdisplay[i] = cnode[3][i];
		histogramweight[i] = length_weight[i];
		if(cnode[3][i] > (INLzonetop + 10.)) totallength_weightfrac += length_weight[i];
		else if(cnode[3][i] >= (INLzonetop - 10.)) INLzonetopfrac += length_weight[i];
		else if(cnode[3][i] < (INLzonebottom - 10.)) totallength_weightfrac += length_weight[i];
		else if(cnode[3][i] <= (INLzonebottom + 10.)) INLzonebottomfrac += length_weight[i];
		else totallength_weightfrac += length_weight[i];
		totallength_weight += length_weight[i];
	}

	if((INLzonetopfrac + INLzonebottomfrac + totallength_weightfrac) - totallength_weight > .001) printf("\n****INLZonelength fractions do not add to totallength!\n***");

	percentINLtop = INLzonetopfrac / totallength_weight * 100;
	percentINLbottom = INLzonebottomfrac / totallength_weight * 100;
	percentneither = totallength_weightfrac / totallength_weight * 100;

	sprintf(fname, "Histogram-Zvessel%03i.out", imain);
	histogram(histogramdisplay, histogramweight, nnod, fname);

	// plot average tissue PO2 vs depth in tissue
	sprintf(fname, "oxy_z_dist%03i.out", imain);
	foxy_z = fopen(fname, "w");
	for (iz = 1; iz <= mzz; iz++) {
		oxy_z[iz] = 0.;
		noxy_z[iz] = 0;
	}
	for (i = 1; i <= nnt; i++) {
		iz = tisspoints[3][i];
		oxy_z[iz] += pt[i][1];
		noxy_z[iz]++;
	}
	for (iz = 1; iz <= mzz; iz++) {
		oxy_z[iz] = oxy_z[iz] / noxy_z[iz];
		zval = (alz / mzz) * (iz - 0.5);
		fprintf(foxy_z, "%f %f\n", zval, oxy_z[iz]);
	}
	fclose(foxy_z);

	free_vector(histogramdisplay, 1, maxdim);
	free_vector(histogramweight, 1, maxdim);
	free_vector(length_weight, 1, nnod);
}