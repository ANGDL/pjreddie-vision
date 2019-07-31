#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // Done
	float s = 0.0;
	for (int i = 0; i != im.w; ++i) {
		for (int j = 0; j != im.h; ++j) {
			for (int k = 0; k != im.c; ++k) {
				s += get_pixel(im, i, j, k);
			}
		}
	}

	for (int i = 0; i != im.w; ++i) {
		for (int j = 0; j != im.h; ++j) {
			for (int k = 0; k != im.c; ++k) {
				float v = get_pixel(im, i, j, k);
				set_pixel(im, i, j, k, v / s);
			}
		}
	}
}

image make_box_filter(int w)
{
    // Done
	image im = make_image(w, w, 1);
	float v = (1.0 / (float)im.w) * (1.0 / (float)im.h);
	for (int i = 0; i != im.w; ++i) {
		for (int j = 0; j != im.h; ++j) {
			for (int k = 0; k != im.c; ++k) {
				set_pixel(im, i, j, k, v);
			}
		}
	}
	//l1_normalize(im);
	return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // Done
	image new_image = make_image(im.w, im.h, im.c);

	int n_filters = filter.c;

	assert(n_filters == im.c || n_filters == 1);
	
	int filter_size = filter.w * filter.h;
	for (int k = 0; k != im.c; ++k) {
		int filter_c = filter.c == 1 ? 0 : k;
		for (int j = 0; j != im.h; ++j) {
			for (int i = 0; i != im.w; ++i) {
				// ¾í»ý
				float v = 0.0;
				for (int m = 0; m != filter.h; ++m) {
					for (int n = 0; n != filter.w; ++n) {
						v += get_pixel(im, i + n - filter.w / 2, j + m -filter.h / 2, k) * get_pixel(filter, n, m, filter_c);
					}
				}
				set_pixel(new_image, i, j, k, v);

			}	
		}
	}

	if(preserve != 0)
		return new_image;

	image avg_img = make_image(new_image.w, new_image.h, 1);
	for (int j = 0; j != new_image.h; ++j) {
		for (int i = 0; i != new_image.w; ++i) {
			float v = 0.0;
			for (int k = 0; k != new_image.c; ++k) {
				v += get_pixel(new_image, i, j, k);
			}

			v /= new_image.c;
			set_pixel(avg_img, i, j, 0, v);
		}
	}
	return avg_img;
}

image make_highpass_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_sharpen_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_emboss_filter()
{
    // TODO
    return make_image(1,1,1);
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    return make_image(1,1,1);
}

image add_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
}

image *sobel_image(image im)
{
    // TODO
    return calloc(2, sizeof(image));
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
