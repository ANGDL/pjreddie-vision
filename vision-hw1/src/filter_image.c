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
			set_pixel(avg_img, i, j, 0, v);
		}
	}
	return avg_img;
}

image make_highpass_filter()
{
    // Done
	float data[3][3] = { {0., -1., 0.},
							 {-1., 4., -1.},
							 {0., -1., 0.} };

	image filter = make_image(3, 3, 1);

	for (int j = 0; j != 3; ++j) {
		for (int i = 0; i != 3; ++i) {
			set_pixel(filter, i, j, 0, data[j][i]);
		}
	}
    return filter;
}

image make_sharpen_filter()
{
    // TODO
	float data[3][3] = { {0., -1., 0.},
							 {-1., 5., -1.},
							 {0., -1., 0.} };

	image filter = make_image(3, 3, 1);

	for (int j = 0; j != 3; ++j) {
		for (int i = 0; i != 3; ++i) {
			set_pixel(filter, i, j, 0, data[j][i]);
		}
	}
	return filter;
}

image make_emboss_filter()
{
	float data[3][3] = { {-2., -1., 0.},
							 {-1., 1., 1.},
							 {0., 1., 2.} };

	image filter = make_image(3, 3, 1);

	for (int j = 0; j != 3; ++j) {
		for (int i = 0; i != 3; ++i) {
			set_pixel(filter, i, j, 0, data[j][i]);
		}
	}
	return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
	/*
	99% of the probability mass for a gaussian is within +/- 3 standard deviations so make the kernel be 6 times the size of sigma. 
	But also we want an odd number, so make it be the next highest odd integer from 6x sigma.
	*/
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

	image filter = make_image(size, size, 1);
	int center = size / 2;

	float s = 0.;
	for (int j = 0; j != filter.h; ++j) {
		for (int i = 0; i != filter.w; ++i) {
			float x = (i - center) * (i - center);
			float y = (j - center) * (j - center);
			float v = exp(-(x + y) / (2 * sigma * sigma)) / (TWOPI * sigma * sigma);
			set_pixel(filter, i, j, 0, v);
			s += v;
		}
	}

	for (int j = 0; j != filter.h; ++j) {
		for (int i = 0; i != filter.w; ++i) {
			float v = get_pixel(filter, i, j, 0);
			set_pixel(filter, i, j, 0, v / s);
		}
	}

    return filter;
}

image add_image(image a, image b)
{
    // Done
	assert(a.w == b.w && a.h == b.h && a.c == b.c);
	image img = make_image(a.w, a.h, a.c);

	for (int k = 0; k != a.c; ++k) {
		for (int j = 0; j != a.h; ++j) {
			for (int i = 0; i != a.w; ++i) {
				float v1 = get_pixel(a, i, j, k);
				float v2 = get_pixel(b, i, j, k);
				set_pixel(img, i, j, k, v1 + v2);
			}
		}
	}

	return img;
}

image sub_image(image a, image b)
{
	// Done
	assert(a.w == b.w && a.h == b.h && a.c == b.c);
	image img = make_image(a.w, a.h, a.c);

	for (int k = 0; k != a.c; ++k) {
		for (int j = 0; j != a.h; ++j) {
			for (int i = 0; i != a.w; ++i) {
				float v1 = get_pixel(a, i, j, k);
				float v2 = get_pixel(b, i, j, k);
				int idx = i + a.w * j + a.w * a.h * k;
				if (idx == 40071) {
					printf("%f,%f\n", v1, v2);
				}
				set_pixel(img, i, j, k, v1 - v2);
			}
		}
	}

	return img;
}

image make_gx_filter()
{
    // TODO
	float data[3][3] = { {-1., 0., 1.},
							 {-2., 0., 2.},
							 {-1., 0., 1.} };

	image filter = make_image(3, 3, 1);

	for (int j = 0; j != 3; ++j) {
		for (int i = 0; i != 3; ++i) {
			set_pixel(filter, i, j, 0, data[j][i]);
		}
	}
	return filter;
}

image make_gy_filter()
{
    // DONE
	float data[3][3] = { {-1., -2., -1.},
								{0., 0., 0.},
								{1., 2., 1.} };

	image filter = make_image(3, 3, 1);

	for (int j = 0; j != 3; ++j) {
		for (int i = 0; i != 3; ++i) {
			set_pixel(filter, i, j, 0, data[j][i]);
		}
	}
	return filter;
}

//min-max normalization
void feature_normalize(image im)
{
	float max_v = -256.0;
	float min_v = 256.0;
	for (int k = 0; k != im.c; ++k) {
		for (int j = 0; j != im.h; ++j) {
			for (int i = 0; i != im.w; ++i) {
				float v = get_pixel(im, i, j, k);
				if (v > max_v)
					max_v = v;
				else if (v < min_v)
					min_v = v;
			}
		}
	}

	for (int i = 0; i != im.w * im.h * im.c; ++i) {
		im.data[i] = (min_v == max_v) ? min_v : (im.data[i] - min_v) / (max_v - min_v);
	}
}

image *sobel_image(image im)
{
    // TODO
	image* img = calloc(2, sizeof(image));
	img[0] = make_image(im.w, im.h, 1);
	img[1] = make_image(im.w, im.h, 1);

	image gx_filter = make_gx_filter();
	image gy_filter = make_gy_filter();

	image gx = convolve_image(im, gx_filter, 0);
	image gy = convolve_image(im, gy_filter, 0);

	save_image(gx, "gx");
	save_image(gy, "gy");

	for (int j = 0; j != im.h; ++j) {
		for (int i = 0; i != im.w; ++i) {
			float v1 = get_pixel(gx, i, j, 0);
			float v2 = get_pixel(gy, i, j, 0);
			float gradient = sqrtf(v1 * v1 + v2 * v2);
			float direction = atan2f(v2, v1);
			set_pixel(img[0], i, j, 0, gradient);
			set_pixel(img[1], i, j, 0, direction);
		}
	}
	
    return img;
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
