#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Comparator for matches
// const void *a, *b: pointers to the matches to compare.
// returns: result of comparison, 0 if same, 1 if a > b, -1 if a < b.
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return  1;
    else return 0;
}

// Helper function to create 2d points.
// float x, y: coordinates of point.
// returns: the point.
point make_point(float x, float y)
{
    point p;
    p.x = x; p.y = y;
    return p;
}

// Place two images side by side on canvas, for drawing matching pixels.
// image a, b: images to place.
// returns: image with both a and b side-by-side.
image both_images(image a, image b)
{
    image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);
    int i,j,k;
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(both, i, j, k, get_pixel(a, i, j, k));
            }
        }
    }
    for(k = 0; k < b.c; ++k){
        for(j = 0; j < b.h; ++j){
            for(i = 0; i < b.w; ++i){
                set_pixel(both, i+a.w, j, k, get_pixel(b, i, j, k));
            }
        }
    }
    return both;
}

// Draws lines between matching pixels in two images.
// image a, b: two images that have matches.
// match *matches: array of matches between a and b.
// int n: number of matches.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i,j;
    for(i = 0; i < n; ++i){
        int bx = matches[i].p.x; 
        int ex = matches[i].q.x; 
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for(j = bx; j < ex + a.w; ++j){
			// 等比三角形
            int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
            set_pixel(both, j, r, 0, i<inliers?0:1);
            set_pixel(both, j, r, 1, i<inliers?1:0);
            set_pixel(both, j, r, 2, 0);
        }
    }
    return both;
}

// Draw the matches with inliers in green between two images.
// image a, b: two images to match.
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// Find corners, match them, and draw them between two images.
// image a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Calculates L1 distance between to floating point arrays.
// float *a, *b: arrays to compare.
// int n: number of values in each array.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(float *a, float *b, int n)
{
    // Done: return the correct number.
	float s = 0.0;
	for (int i = 0; i != n; ++i) {
		float v = a[i] - b[i];
		s += (v > 0) ? v : -v;
	}
    return s;
}

// Finds best matches between descriptors of two images.
// descriptor *a, *b: array of descriptors for pixels in two images.
// int an, bn: number of descriptors in arrays a and b.
// int *mn: pointer to number of matches found, to be filled in by function.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i,j;

    // We will have at most an matches.
    *mn = an;
    match *m = (match*)calloc(an, sizeof(match));
    for(j = 0; j < an; ++j){
        // TODO: for every descriptor in a, find best match in b.
        // record ai as the index in *a and bi as the index in *b.
        int bind = 0; // <- find the best match
		float distance = 99999.;
		for (i = 0; i < bn; ++i) {
			assert(a[j].n == b[i].n);
			float d = l1_distance(a[j].data, b[i].data, a[j].n);
			if (d < distance) {
				distance = d;
				bind = i;
			}
		}
        m[j].ai = j;
        m[j].bi = bind; // <- should be index in b.
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = distance; // <- should be the smallest L1 distance!
    }

    int count = 0;
    int *seen = calloc(bn, sizeof(int));
    // DONE: we want matches to be injective (one-to-one).
    // Sort matches based on distance using match_compare and qsort.
    // Then throw out matches to the same element in b. Use seen to keep track.
    // Each point should only be a part of one match.
    // Some points will not be in a match.
    // In practice just bring good matches to front of list, set *mn.

	qsort(m, an, sizeof(match), match_compare);

	//for (int t = 0; t < an; ++t) {
	//	printf("%d,", m[t].bi);
	//}

	i = 0, j = i + 1;
	if (seen) {
		++seen[m[i].bi];

		while (j < an)
		{
			if (seen[m[j].bi] == 0) {
				if (i + 1 != j) {
					m[i + 1] = m[j];
				}
				++seen[m[j].bi];
				++i;
				++j;
				++count;
			}
			else if (seen[m[j].bi] > 0) {
				++j;
			}
		}

		*mn = count + 1;
		free(seen);
	}

    return m;
}

// Apply a projective transformation to a point.
// matrix H: homography to project point.
// point p: point to project.
// returns: point projected using the homography.
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    // DONE: project point p with homography H.
    // Remember that homogeneous coordinates are equivalent up to scalar.
    // Have to divide by.... something...
    point q = make_point(0, 0);
	c.data[0][0] = p.x;
	c.data[1][0] = p.y;
	c.data[2][0] = 1;

	matrix r = matrix_mult_matrix(H, c);

	q.x = (float) (r.data[0][0] / r.data[2][0]);
	q.y = (float) (r.data[1][0] / r.data[2][0]);

	free_matrix(c);
	free_matrix(r);

    return q;
}

// Calculate L2 distance between two points.
// point p, q: points.
// returns: L2 distance between them.
float point_distance(point p, point q)
{
    // DONE: should be a quick one.
	float a = p.x - q.x;
	float b = p.y - q.y;
	float d = sqrtf(a * a + b * b);
    return d;
}

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// matrix H: homography between coordinate systems.
// match *m: matches to compute inlier/outlier.
// int n: number of matches in m.
// float thresh: threshold to be an inlier.
// returns: number of inliers whose projected point falls within thresh of
//          their match in the other image. Should also rearrange matches
//          so that the inliers are first in the array. For drawing.
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i = 0;
    int count = 0;
	if (H.data == NULL || m == NULL) {
		return count;
	}
    // DONE: count number of matches that are inliers
    // i.e. distance(H*p, q) < thresh
    // Also, sort the matches m so the inliers are the first 'count' elements.
	
	
	//printf("model_inliers distance:\n---------------------\n");
	//while (i < n)
	//{
	//	point p = project_point(H, m[i].p);
	//	point q = m[i].q;
	//	float d = point_distance(p, q);
	//	printf("%f, ", d);
	//	++i;
	//}
	//printf("\n---------------------\n");

	i = 0;
	int j = i + 1 ;
	while (j <= n)
	{
		point p = project_point(H, m[i].p);
		point q = m[i].q;
		float d = point_distance(p, q);
		if (d >= thresh) {
			if (j == n) {
				break;
			}
			match t = m[j];		
			m[j] = m[i];
			m[i] = t;
			++j;
		}
		else {
			++i;
			++j;
			++count;
		}
	}

	//printf("\n");

	//for (i = 0; i < count; i++)
	//{
	//	point p = project_point(H, m[i].p);
	//	point q = m[i].q;
	//	float d = point_distance(p, q);
	//	printf("%f, ", d);
	//}
	//printf("\n---------------------\n");

    return count;
}

// Randomly shuffle matches for RANSAC.
// match *m: matches to shuffle in place.
// int n: number of elements in matches.
void randomize_matches(match *m, int n)
{
    // DONE: implement Fisher-Yates to shuffle the array.

	for (int i = n - 1; i > 0; i--) {
		int j = rand() % (i + 1);
		match swap = m[i];
		m[i] = m[j];
		m[j] = swap;
	}
}

// Computes homography between two images given matching pixels.
// match *matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
matrix compute_homography(match *matches, int n)
{
    matrix M = make_matrix(n*2, 8);
    matrix b = make_matrix(n*2, 1);

	assert(n > 0 && matches != NULL);
	//if (n <= 0 || matches == NULL) {
	//	printf("compute_homography input error for n=%d, matches=%s", n, matches);
	//	return make_matrix(3, 3);
	//}

    int i;
    for(i = 0; i < n; ++i){
        double x  = matches[i].p.x;
        double xp = matches[i].q.x;
        double y  = matches[i].p.y;
        double yp = matches[i].q.y;
        // TODO: fill in the matrices M and b.
		int r1 = i * 2;
		int r2 = r1 + 1;

		M.data[r1][0] = x;
		M.data[r1][1] = y;
		M.data[r1][2] = 1;
		M.data[r1][6] = -x * xp;
		M.data[r1][7] = -y * xp;
		b.data[r1][0] = xp;

		M.data[r2][3] = x;
		M.data[r2][4] = y;
		M.data[r2][5] = 1;
		M.data[r2][6] = -x * yp;
		M.data[r2][7] = -y * yp;
		b.data[r2][0] = yp;
    }
    matrix a = solve_system(M, b);
    free_matrix(M); free_matrix(b); 

    // If a solution can't be found, return empty matrix;
    matrix none = {0};
    if(!a.data) return none;

    matrix H = make_matrix(3, 3);
    // DONE: fill in the homography H based on the result in a.
	for (i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			int idx = i * 3 + j;
			if (idx == 8) {
				H.data[i][j] = 1;
			}
			else {
				H.data[i][j] = a.data[idx][0];
			}
		}
	}

	// TEST
	//printf("H:\n");
	//for (i = 0; i < 3; ++i) {
	//	for (int j = 0; j < 3; ++j) {
	//		int idx = i * 3 + j;
	//		printf("%f ", H.data[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("---------------------\n");
	//for (i = 0; i < 3; ++i) {
	//	for (int j = 0; j < 3; ++j) {
	//		int idx = i * 3 + j;
	//		if (idx < 8)
	//		{
	//			printf("%f ", a.data[idx][0]);
	//		}
	//	}
	//	printf("\n");
	//}
	//printf("---------------------\n");

    free_matrix(a);
    return H;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// match *m: set of matches.
// int n: number of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
matrix RANSAC(match* m, int n, float thresh, int k, int cutoff)
{
	int best = 0;
	matrix Hb = make_translation_homography(256, 0);
	// TODO: fill in RANSAC algorithm.
	// for k iterations:
	//     shuffle the matches
	//     compute a homography with a few matches (how many??)
	//     if new homography is better than old (how can you tell?):
	//         compute updated homography using all inliers
	//         remember it and how good it is
	//         if it's better than the cutoff:
	//             return it immediately
	// if we get to the end return the best homography
	for (int i = 0; i < k; i++) {
		randomize_matches(m, n);
		matrix H = compute_homography(m, 4);
		int curr_inliers = model_inliers(H, m, n, thresh);
		if (curr_inliers > best) {
			best = curr_inliers;
			free_matrix(Hb);
			Hb = H;
			if (best > cutoff) {
				return Hb;
			}
		}
	}
	return Hb;
}

// Stitches two images together using a projective transformation.
// image a, b: images to stitch.
// matrix H: homography from image a coordinates to image b coordinates.
// returns: combined image stitched together.
image combine_images(image a, image b, matrix H)
{
    matrix Hinv = matrix_invert(H);

    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // Can disable this if you are making very big panoramas.
    // Usually this means there was an error in calculating H.
    if(w > 7000 || h > 7000){
        fprintf(stderr, "output too big, stopping\n");
        return copy_image(a);
    }

    int i,j,k;
    image c = make_image(w, h, a.c);
    
    // Paste image a into the new image offset by dx and dy.
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                // DONE: fill in.
				int nx = i - dx;
				int ny = j - dy;
				float v = get_pixel(a, i, j, k);
				set_pixel(c, nx, ny, k, v);
            }
        }
    }

    // DONE: Paste in image b as well.
    // You should loop over some points in the new image (which? all?)
    // and see if their projection from a coordinates to b coordinates falls
    // inside of the bounds of image b. If so, use bilinear interpolation to
    // estimate the value of b at that projection, then fill in image c.
	for (k = 0; k < c.c; ++k) {
		for (j = 0; j < c.h; ++j) {
			for (i = 0; i < c.w; ++i) {
				// find the point in b projected using the homography
				point b_point = project_point(H, make_point((float)i, (float)j));
				// if projected point in b.w and b.h
				if (b_point.x >=0 && b_point.x < b.w && b_point.y >= 0 && b_point.y < b.h) {
					// calculate the pixel
					float b_pix = bilinear_interpolate(b, b_point.x, b_point.y, k);
					// 计算偏移后的坐标
					int nx = i - dx;
					int ny = j - dy;
					if(nx >= 0 && nx < c.w && ny >= 0 && ny < c.h)
						set_pixel(c, nx, ny, k, b_pix);
				}
			}
		}
	}
    return c;
}

// Create a panoramam between two images.
// image a, b: images to stitch together.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
// float inlier_thresh: threshold for RANSAC inliers. Typical: 2-5
// int iters: number of RANSAC iterations. Typical: 1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;
    
    // Calculate corners and descriptors
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    // Find matches
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    // Run RANSAC to find the homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    if(0){
        // Mark corners and matches between images
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);

    // Stitch the images together with the homography
    image comb = combine_images(a, b, H);
    return comb;
}

// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
image cylindrical_project(image im, float f)
{
    //DONE: project image onto a cylinder
    image c = copy_image(im);
	int xc = im.w / 2;
	int yc = im.h / 2;

	for (int k = 0; k != im.c; ++k) {
		for (int j = 0; j != im.h; j++)
		{
			for (int i = 0; i != im.w; ++i)
			{
				float theta = (i - xc) / f;
				float h = (j - yc) / f;
				float x = sinf(theta);
				float y = h;
				float z = cosf(theta);
				x = f * x / z + xc;
				y = f * y / z + yc;
				float v = get_pixel(im, i, j, k);
				set_pixel(c, (int)(x), (int)(y), k, v);
			}
		}
	}

    return c;
}
