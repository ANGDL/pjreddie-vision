#include <math.h>
#include <string.h>
#include "image.h"
#include "test.h"
#include "args.h"

void easy_panorama() {

	image a = load_image("data/Rainier1.png");
	image b = load_image("data/Rainier2.png");
	image m = panorama_image(a, b, 2, 50, 3, 5, 1000, 50);
	save_image(m, "matches");
}

void rainier_panorama() {
	image im1 = load_image("data/Rainier1.png");
	image im2 = load_image("data/Rainier2.png");
	image im3 = load_image("data/Rainier3.png");
	image im4 = load_image("data/Rainier4.png");
	image im5 = load_image("data/Rainier5.png");
	image im6 = load_image("data/Rainier6.png");
	image pan = panorama_image(im1, im2, 2, 50, 3, 5, 1000, 50);
	save_image(pan, "rainier_panorama_1");
	image pan2 = panorama_image(pan, im5, 2, 50, 3, 5, 1000, 50);
	save_image(pan2, "rainier_panorama_2");
	image pan3 = panorama_image(pan2, im6, 2, 50, 3, 5, 1000, 50);
	save_image(pan3, "rainier_panorama_3");
	image pan4 = panorama_image(pan3, im3, 2, 50, 3, 5, 1000, 50);
	save_image(pan4, "rainier_panorama_4");
	image pan5 = panorama_image(pan4, im4, 2, 50, 3, 5, 1000, 50);
	save_image(pan5, "rainier_panorama_5");
}

int main(int argc, char **argv)
{
    //char *in = find_char_arg(argc, argv, "-i", "data/dog.jpg");
    //char *out = find_char_arg(argc, argv, "-o", "out");
    ////float scale = find_float_arg(argc, argv, "-s", 1);
    //if(argc < 2){
    //    printf("usage: %s [test | grayscale]\n", argv[0]);  
    //} else if (0 == strcmp(argv[1], "test")){
    //    run_tests();
    //} else if (0 == strcmp(argv[1], "grayscale")){
    //    image im = load_image(in);
    //    image g = rgb_to_grayscale(im);
    //    save_image(g, out);
    //    free_image(im);
    //    free_image(g);
    //}

	rainier_panorama();
	return 0;
}
