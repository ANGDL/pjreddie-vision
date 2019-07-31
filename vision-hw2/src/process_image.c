#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // Done Fill this in
	// bounds checking
	if (x < 0) x = 0;
	if (x >= im.w) x = im.w - 1;
	if (y < 0) y = 0;
	if (y >= im.h) y = im.h - 1;
	if (c < 0) c = 0;
	if (c >= im.c) c = im.c - 1;

	// get the index
	int idx = x + im.w * y + im.w * im.h * c;
	assert(idx < im.w * im.h * im.c);
    return im.data[idx];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // Done Fill this in
	// get the index
	int idx = x + im.w * y + im.w * im.h * c;
	assert(idx < im.w * im.h * im.c);
	im.data[idx] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // Done Fill this in
	memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // Done
	for (int j = 0; j != im.h; ++j) {
		for (int i = 0; i != im.w; ++i) {
			float v = 0.299 * get_pixel(im, i, j, 0) + 0.587 * get_pixel(im, i, j, 1) + 0.114 * get_pixel(im, i, j, 2);
			set_pixel(gray, i, j, 0, v);
		}
	}
    return gray;
}

void shift_image(image im, int c, float v)
{
    // Done
	for (int j = 0; j != im.h; ++j) {
		for (int i = 0; i != im.w; ++i) {
			float new_v = get_pixel(im, i, j, c) + v;
			set_pixel(im, i, j, c, new_v);
		}
	}
}

void clamp_image(image im)
{
	// Done
	int size_column = im.w;
	int size_row = im.h;
	int size_channel = im.c;

	for (int k = 0; k < size_channel; ++k) {
		for (int j = 0; j < size_row; ++j) {
			for (int i = 0; i < size_column; ++i) {
				int index = i + size_column * j + size_column * size_row * k;
				float value = im.data[index];
				if (value > 1) im.data[index] = 1;
				else if (value < 0) im.data[index] = 0;
			}
		}
	}
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
	assert(3 == im.c);
    // Done
	for (int j = 0; j != im.h; ++j) {
		for (int i = 0; i != im.w; ++i) {
			float r = get_pixel(im, i, j, 0);
			float g = get_pixel(im, i, j, 1);
			float b = get_pixel(im, i, j, 2);

			float v = three_way_max(r, g, b);
			float m = three_way_min(r, g, b);

			float c = v - m;
			float s = (v > 0.0)? c / v : 0.0;

			float h = 0.0; // c == 0

			if (c > 0.0) {
				if (v == r) {
					h = (g - b) / c;
				}
				else if (v == g) {
					h = (b - r) / c + 2.;
				}
				// v == b
				else if (v == b) {
					h = (r - g) / c + 4.;
				}
			}

			if (h < 0) {
				h = (h / 6. + 1.);
			}
			else {
				h /= 6.;
			}

			set_pixel(im, i, j, 0, h);
			set_pixel(im, i, j, 1, s);
			set_pixel(im, i, j, 2, v);
		}
	}
}

void hsv_to_rgb(image im)
{
	assert(3 == im.c);
	// Done
	for (int y = 0; y < im.h; y++)
	{
		for (int x = 0; x < im.w; x++)
		{
			float h = get_pixel(im, x, y, 0);
			float s = get_pixel(im, x, y, 1);
			float v = get_pixel(im, x, y, 2);

			float r = 0;
			float g = 0;
			float b = 0;

			int i = (int)(h * 6);
			float f = (h * 6) - (float)i;
			float p = v * (1 - s);
			float q = v * (1 - s * f);
			float t = v * (1 - s * (1 - f));

			if (i == 0)
			{
				r = v;
				g = t;
				b = p;
			}
			else if (i == 1)
			{
				r = q;
				g = v;
				b = p;
			}
			else if (i == 2)
			{
				r = p;
				g = v;
				b = t;
			}
			else if (i == 3)
			{
				r = p;
				g = q;
				b = v;
			}
			else if (i == 4)
			{
				r = t;
				g = p;
				b = v;
			}
			else if (i == 5)
			{
				r = v;
				g = p;
				b = q;
			}

			set_pixel(im, x, y, 0, r);
			set_pixel(im, x, y, 1, g);
			set_pixel(im, x, y, 2, b);
		}
	}
}
