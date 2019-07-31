#include <math.h>
#include "image.h"
#include<assert.h>

float nn_interpolate(image im, float x, float y, int c)
{
    // Done
	int i = rint(x);
	int j = rint(y);
	float v = get_pixel(im, i, j, c);
    return v;
}

image nn_resize(image im, int w, int h)
{
	image new_image = make_image(w, h, im.c);

	float width_ratio = ((float)im.w) / w;
	float height_ratio = ((float)im.h) / h;

	for (int y = 0; y != new_image.h; ++y)
	{
		for (int x = 0; x != new_image.w; ++x)
		{
			for (int z = 0; z != new_image.c; ++z)
			{
				// 源图像和目标图像几何中心的对齐 
				// SrcX = (dstX + 0.5) * (srcWidth / dstWidth) - 0.5
				//	SrcY = (dstY + 0.5) * (srcHeight / dstHeight) - 0.5

				float scaled_x = x * width_ratio + (width_ratio * 0.5 - 0.5);
				float scaled_y = y * height_ratio + (height_ratio * 0.5 - 0.5);
				float value = nn_interpolate(im, scaled_x, scaled_y, z);
				set_pixel(new_image, x, y, z, value);
			}
		}
	}

	return new_image;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO

	float d1 = x - floorf(x);
	float d2 = ceil(x) - x;

	float d3 = y - floorf(y);
	float d4 = ceil(y) - y;

	int x0 = floorf(x);
	int y0 = floorf(y);
	int x1 = ceil(x);
	int y1 = ceil(y);

	float q1 = d4 * get_pixel(im, x0, y0, c) + d3 * get_pixel(im, x0, y1, c);
	float q2 = d4 * get_pixel(im, x1, y0, c) + d3 * get_pixel(im, x1, y1, c);

	float v = d2 * q1 + d1 * q2;

    return v;
}

image bilinear_resize(image im, int w, int h)
{
	image new_image = make_image(w, h, im.c);

	float width_ratio = ((float)im.w) / w;
	float height_ratio = ((float)im.h) / h;

	for (int k = 0; k != new_image.c; ++k) {
		for (int j = 0; j != h; ++j) {
			for (int i = 0; i != w; ++i) {
				// 源图像和目标图像几何中心的对齐 
				// SrcX = (dstX + 0.5) * (srcWidth / dstWidth) - 0.5
				//	SrcY = (dstY + 0.5) * (srcHeight / dstHeight) - 0.5

				float scaled_x = i * width_ratio + (width_ratio * 0.5 - 0.5);
				float scaled_y = j * height_ratio + (height_ratio * 0.5 - 0.5);
				float value = bilinear_interpolate(im, scaled_x, scaled_y, k);
				set_pixel(new_image, i, j, k, value);
			}
		}
	}

	return new_image;
}

