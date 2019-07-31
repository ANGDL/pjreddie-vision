#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
	int size = sigma * 6;
	while (1)
	{
		if (size % 2 == 1) {
			break;
		}
		else {
			++size;
		}
	}

	image filter = make_image(size, 1, 1);
	int center = size / 2;

	float s = 0.;

	for (int i = 0; i != filter.w; ++i) {
		float x = (i - center) * (i - center);
		float v = exp(-x / (2 * sigma * sigma)) / (TWOPI * sigma * sigma);
		set_pixel(filter, i, 0, 0, v);
		s += v;
	}

	for (int i = 0; i != filter.w; ++i) {
		float v = get_pixel(filter, i, 0, 0);
		set_pixel(filter, i, 0, 0, v / s);
	}

	return filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // Done: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
		image g = make_1d_gaussian(sigma);
		image s = convolve_image(im, g, 1);
		free_image(g);
		return s;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
	// 1. calculating derivatives
	if (3 == im.c)
		im = rgb_to_grayscale(im);

	int img_size = im.w * im.h;

	image gx_filter = make_gx_filter();
	image gy_filter = make_gy_filter();

	image gx_img = convolve_image(im, gx_filter, 0);
	image gy_img = convolve_image(im, gy_filter, 0);

	image ix_img = make_image(im.w, im.h, 1);
	image iy_img = make_image(im.w, im.h, 1);
	image ixy_img = make_image(im.w, im.h, 1);

	for (int j = 0; j != im.h; ++j) {
		for (int i = 0; i != im.w; ++i) {
			float ix = get_pixel(gx_img, i, j, 0) /** get_pixel(im, i, j, 0)*/;
			float iy = get_pixel(gy_img, i, j, 0) /** get_pixel(im, i, j, 0)*/;
			float ixy = ix * iy;
			set_pixel(ix_img, i, j, 0, ix * ix);
			set_pixel(iy_img, i, j, 0, iy * iy);
			set_pixel(ixy_img, i, j, 0, ixy);
		}
	}

	//2. corresponding measures
	image sx_img = smooth_image(ix_img, sigma);
	image sy_img = smooth_image(iy_img, sigma);
	image sxy_img = smooth_image(ixy_img, sigma);

	float* p = S.data;
	memmove(p, sx_img.data, img_size * sizeof(float));
	p += img_size;
	memmove(p, sy_img.data, img_size * sizeof(float));
	p += img_size;
	memmove(p, sxy_img.data, img_size * sizeof(float));

	free_image(gx_filter);
	free_image(gy_filter);
	free_image(gx_img);
	free_image(gy_img);
	free_image(ix_img);
	free_image(iy_img);
	free_image(ixy_img);
	
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
	for (int j = 0; j != R.h; ++j) {
		for (int i = 0; i != R.w; ++i) {
			float sx = get_pixel(S, i, j, 0);
			float sy = get_pixel(S, i, j, 1);
			float sxy = get_pixel(S, i, j, 2);

			float det = sx * sy - sxy * sxy;
			float trace = sx + sy;

			float h = det - 0.06 * trace * trace;
			set_pixel(R, i, j, 0, h);
		}
	}
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
	assert(im.c == 1);
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
	for (int j = 0; j != r.h; ++j) {
		for (int i = 0; i != r.w; ++i) {
			float v = get_pixel(r, i, j, 0);
			for (int m = 0; m != w; ++m) {
				for (int n = 0; n != w; ++n) {
					float nv = get_pixel(r, i + m - 2 / w, j + n - 2 / w, 0);
					if (nv > v) {
						set_pixel(r, i, j, 0, -99999.);
						break;
					}
				}
			}
		}
	}

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0; // change this
	for (int i = 0; i != Rnms.w * Rnms.h * Rnms.c; ++i) {
		float v = Rnms.data[i];
		if (v > thresh) {
			++count;
		}
	}
    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
	count = 0;
	for (int i = 0; i != Rnms.w * Rnms.h * Rnms.c; ++i) {
		float v = Rnms.data[i];
		if (v > thresh) {
			d[count++] = describe_index(im, i);
		}
	}

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
