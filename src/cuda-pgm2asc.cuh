#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "unicode_defs.h"
#include "pgm2asc.h"
#include "pnm.h"
#include "gocr.h"


typedef struct box_d box_d;

struct box_d {
	int x0,x1,y0,y1;
	int m1,m2,m3,m4;
	int wac_0, num_ac;
	wchar_t c;
	struct box *box;
};

struct job_d {
	struct {
		int n_run;
	} tmp;

	struct {
		int cs;
		int certainty;
	} cfg;
};

struct return_element {
	box_d box2;
	box_d box4;
	wchar_t bc;
	int j_max;
	int dist;
};


struct return_element *deviceFuncCall (job_t *job, pix *pp);