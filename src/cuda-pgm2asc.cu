extern "C" {
    #include "cuda-pgm2asc.cuh"
};


#define TREE_ARRAY_SIZE 1024
#define Nfilt3 6


void rec_generate_tree_d (char *tree, const char *filter, int i, int n)
{
    assert (i >= 0 && i <= 9);
    assert (n < TREE_ARRAY_SIZE);
	
	if (i == 9) {
		if (filter[4] == 0)
			tree[n] = 2;
		else
			tree[n] = 1;
		
		return;
	}
	
	if (n != -1)
		tree[n] = 1;
	if (filter[i] == 0)
		rec_generate_tree_d (tree, filter, i + 1, n * 2 + 2);
	else if (filter[i] == 1)
		rec_generate_tree_d (tree, filter, i + 1, n * 2 + 3);
	else {
		rec_generate_tree_d (tree, filter, i + 1, n * 2 + 2);
		rec_generate_tree_d (tree, filter, i + 1, n * 2 + 3);
	}
}


__device__
int getpixel_d (pix *p, int x, int y, struct job_d *job, char *tree)
{
	if ( x < 0 || y < 0 || x >= p->x || y >= p->y ) 
	  return 255 & ~7;
  
	if (job->tmp.n_run > 0) {
		int pixel_val = pixel_atp (p, x, y) & ~7;
		int n = -1;
	
	#define IS_BLACK(_dx,_dy) !(pixel_atp (p, x + (_dx), y + (_dy)) >> 7)
	#define IS_WHITE(_dx,_dy) (pixel_atp (p, x + (_dx), y + (_dy)) >> 7)
	#define GO_LEFT n = n * 2 + 2
	#define GO_RIGHT n = n * 2 + 3
	#define CHECK_NO_MATCH if (tree[n] == 0) return pixel_val
	
		if (y == 0) {
			n = 13;
		} else {
			if (x == 0 || IS_BLACK (-1, -1)) 
				GO_RIGHT;
			else  
				GO_LEFT;
	
			if (IS_WHITE (0, -1)) 
				GO_LEFT;
			else  
				GO_RIGHT;
			CHECK_NO_MATCH;
	
			if (x + 1 == p->x || IS_BLACK (+1, -1))
				GO_RIGHT;
			else 
				GO_LEFT;
			CHECK_NO_MATCH;
		}
	
	
		if (x == 0 || IS_BLACK (-1, 0)) 
			GO_RIGHT;
		else 
			GO_LEFT;
		CHECK_NO_MATCH;
	
	
		if (IS_WHITE (0, 0))
			GO_LEFT;
		else
			GO_RIGHT;
		CHECK_NO_MATCH;
	
		if (x + 1 == p->x || IS_BLACK (+1, 0)) 
			GO_RIGHT;
		else 
			GO_LEFT;
		CHECK_NO_MATCH;
	
		if (y + 1 == p->y) {
			n = 8 * n + 21;
		} else {
			if (x == 0 || IS_BLACK (-1, +1)) 
				GO_RIGHT;
			else 
				GO_LEFT;
			CHECK_NO_MATCH;
	
			if (IS_WHITE (0, 1)) 
				GO_LEFT;
			else  
				GO_RIGHT;
			CHECK_NO_MATCH;
	
			if (x + 1 == p->x || IS_BLACK (+1, +1)) 
				GO_RIGHT;
			else 
				GO_LEFT;
		}
	
		CHECK_NO_MATCH;
	
		if (tree[n] == 1) {
			return job->cfg.cs;
		} else {
			return 0;
		}

	}

	return (pixel_atp (p,x,y) & ~7);
}


__device__
int distance_d( pix *p1, box_d *box1, pix *p2, box_d *box2, struct job_d *job, char *tree)
{   
	int rc=0,x,y,v1,v2,i1,i2,rgood=0,rbad=0,x1,y1,x2,y2,dx,dy,dx1,dy1,dx2,dy2;
	int cs = job->cfg.cs;

	x1=box1->x0;
	y1=box1->y0;
	x2=box2->x0;
	y2=box2->y0;
  
	dx1=box1->x1 - box1->x0 + 1;
	dx2=box2->x1 - box2->x0 + 1;
	dx=((dx1>dx2)?dx1:dx2);
  
	dy1=box1->y1 - box1->y0 + 1;
	dy2=box2->y1 - box2->y0 + 1;
	dy=((dy1>dy2)?dy1:dy2);
  
	if(abs(dx1-dx2)>1+dx/16 || abs(dy1-dy2)>1+dy/16) return 100;
	// compare relations to baseline and upper line
	if(2*box1->y1>box1->m3+box1->m4 && 2*box2->y1<box2->m3+box2->m4) rbad += 128;
	if(2*box1->y0>box1->m1+box1->m2 && 2*box2->y0<box2->m1+box2->m2) rbad += 128;
	// compare pixels
	for(y = 0; y < dy; y++) {
	  for(x = 0; x < dx; x++) {	// try global shift too ???
		v1 = ((getpixel_d (p1, x1 + x, y1 + y, job, tree) < cs)?1:0);
		i1=8;	// better gray?
		
		v2 = ((getpixel_d (p2, x2 + x, y2 + y, job, tree) < cs)?1:0);
		i2=8;	// better gray?
		
		if(v1 == v2) {
		  rgood += 8;
		  continue;
		} // all things are right!
		// what about different pixel???
		// test overlap of 8 surounding pixels ??? bad if two nb. are bad
		v1=-1;

		for(i1=-1; i1 < 2; i1++) {
		  for(i2=-1; i2 < 2; i2++) {
			if(i1!=0 || i2!=0){
			  if( ((getpixel_d(p1,x1+x+i1*(1+dx/32),y1+y+i2*(1+dy/32), job, tree)<cs)?1:0)
			  !=((getpixel_d(p2,x2+x+i1*(1+dx/32),y2+y+i2*(1+dy/32), job, tree)<cs)?1:0) ) v1++;
			}
		  }
		}

		if (v1 > 0) rbad+=16*v1;
		else rbad++;    
	  }
	}
  
	if(rgood + rbad) rc = (100*rbad+(rgood+rbad-1))/(rgood+rbad);
	else rc = 99;
	
	return rc;
  }


__global__
void deviceFunc (int n, box_d *boxArr, struct job_d *job, pix *pp, struct return_element *returnArr, char *tree)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	if (i < n && j < n) {
		struct box_d *box3 = boxArr + j;
		struct box_d *box2 = boxArr + i;
		int wac = ((box3->num_ac > 0)?box3->wac_0:100);
		int d;

		if (box2->c == UNKNOWN || (box2->num_ac > 0 && box2->wac_0 < 97)) {
			if (box2->y1 - box2->y0 > 4 && box2->x1 - box2->x0 > 1) {

				int *dist = &(returnArr[i].dist);
				int *j_max = &(returnArr[i].j_max);

				if (box3 == box2 || box3->c == UNKNOWN || wac < job->cfg.certainty);
				else if (box2->y1 - box2->y0 < 5 || box2->x1 - box2->x0 < 3);
				else {

					d = distance_d (pp, box2, pp, box3, job, tree);

					atomicMin (dist, d);

					__syncthreads ();

					if (d == *dist) {
						atomicMax (j_max, j);
					}

					__syncthreads ();

					if (j == *j_max) {
						returnArr[i].box2 = *box2;
						returnArr[i].box4 = *box3;
						returnArr[i].bc = box3->c;
					}
				}
			}
		}
	}
}

struct return_element *deviceFuncCall (job_t *job, pix *pp)
{
	const char filt3[Nfilt3][9] = { 
		{0,0,0, 0,0,1, 1,0,0},
		{0,0,0, 1,0,1, 0,0,0},
		{1,0,0, 0,0,1, 0,0,0},
		{1,1,0, 0,1,0, 2,1,1},
		{0,0,1, 0,0,0, 2,1,0},
		{0,1,0, 0,0,0, 1,2,0}
	};
	
	char tree[TREE_ARRAY_SIZE];

	memset (tree, 0, sizeof(tree));

	for (int f = 0; f < Nfilt3; f++) {
		const char * filter = filt3[f];
		rec_generate_tree_d (tree, filter, 0, -1);
	}

	struct timeval stop, start;

	int n = job->res.boxlist.n;

	box_d *boxArr = (box_d *) malloc (n * sizeof (box_d)); 
	struct job_d *jobArr = (struct job_d *) malloc (sizeof (struct job_d));
	pix *pixArr = (pix *) malloc (sizeof (pix));
	struct return_element *returnArr = (struct return_element *) malloc (n * sizeof (struct return_element));

	int i = 0;
	for_each_data (&(job->res.boxlist)) {
		struct box *tmpBox = (struct box *) list_get_current (&(job->res.boxlist));

		boxArr[i] = {
			tmpBox->x0, tmpBox->x1, tmpBox->y0, tmpBox->y1,
			tmpBox->m1, tmpBox->m2, tmpBox->m3, tmpBox->m4,
			tmpBox->wac[0], tmpBox->num_ac, tmpBox->c,
			tmpBox
		};

		i++;
	} end_for_each (&(job->res.boxlist));

	*jobArr = {
		{
			job->tmp.n_run
		},
		{
			job->cfg.cs,
			job->cfg.certainty
		}
	};

	unsigned char *p_d;
	cudaMalloc (&p_d, pp->x * pp->y);
	cudaMemcpy (p_d, pp->p, pp->x * pp->y, cudaMemcpyHostToDevice);
	*pixArr = {
		p_d, pp->x, pp->y, pp->bpp
	};

	for (int i = 0; i < n; i++) {
		returnArr[i].j_max = -1;
		returnArr[i].dist = 1000;
	}

	box_d *boxArr_d;
	struct job_d *jobArr_d;
	pix *pixArr_d;
	struct return_element *returnArr_d;
	char *treeArr_d;

	cudaMalloc (&boxArr_d, n * sizeof (box_d));
	cudaMalloc (&jobArr_d, sizeof (struct job_d));
	cudaMalloc (&pixArr_d, sizeof (pix));
	cudaMalloc (&returnArr_d, n * sizeof (struct return_element));
	cudaMalloc (&treeArr_d, sizeof (tree));

	cudaMemcpy (boxArr_d, boxArr, n * sizeof (box_d), cudaMemcpyHostToDevice);
	cudaMemcpy (jobArr_d, jobArr, sizeof (struct job_d), cudaMemcpyHostToDevice);
	cudaMemcpy (pixArr_d, pixArr, sizeof (pix), cudaMemcpyHostToDevice);
	cudaMemcpy (treeArr_d, tree, sizeof (tree), cudaMemcpyHostToDevice);
	cudaMemcpy (returnArr_d, returnArr, n * sizeof (struct return_element), cudaMemcpyHostToDevice);


    dim3 threadsPerBlock (16, 16);
    dim3 numBlocksTemp (0, 0);

    if ((n / threadsPerBlock.x) * threadsPerBlock.x < n) {
        numBlocksTemp.x = n / threadsPerBlock.x + 1;
    }
    if ((n / threadsPerBlock.y) * threadsPerBlock.y < n) {
        numBlocksTemp.y = n / threadsPerBlock.y + 1;
    }

    dim3 numBlocks (numBlocksTemp.x, numBlocksTemp.y);


	fprintf (stderr, "\n==========\ngpu compute\n");
	fprintf (stderr, "\ngridDim %d, %d\nblockDim %d, %d\n", numBlocks.x, numBlocks.y, threadsPerBlock.y, threadsPerBlock.y);
	gettimeofday (&start, NULL);
	
	deviceFunc <<<numBlocks, threadsPerBlock>>> (n, boxArr_d, jobArr_d, pixArr_d, returnArr_d, treeArr_d);
	cudaDeviceSynchronize ();

	gettimeofday (&stop, NULL);
    fprintf (stderr, "\ntook %lu us\n==========\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);


	fprintf (stderr, "\n\n==========\ndevice to host memcpy\n");
	gettimeofday (&start, NULL);

	cudaMemcpy (returnArr, returnArr_d, n * sizeof (struct return_element), cudaMemcpyDeviceToHost);

	gettimeofday (&stop, NULL);
    fprintf (stderr, "\ntook %lu us\n==========\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);

	return returnArr;
}