#pragma once

#include <map>
#include <Image.hpp>
#include <Math.hpp>
#include <Area.hpp>
#include <Point.hpp>


class ImageUtil{
private:
	float static laplacian_filter[9];
	int static min_calc(int *p);
	int static Lookup_update(int i,int * label);
	void static filtering(Image &image, float *filter, int filter_size[2], int direction[2]);
	float static convolution(Image &image, float *filter, int x, int y, int filter_size[2], int direction[2]);
	double static derivativesGaussian(int x, int y, float sigma);
	double static Gaussian(int x, int y, float sigma);
public:
	void static Shrink(Image &image);
	void static Expand(Image &image);
	void static Opening(Image &image, int level);
	void static Closing(Image &image, int level);
	void static TopHat(Image &image,int level);
	void static BlackHat(Image &image,int level);
	int static Thinning(Image &image);
	void static UnSharpMasking(Image &image,int level);
	void static LowResolution(Image &image,int level);
	void static colorRegionSplit(Image &image, int threshold = 100);
	void static getHistogram(Image &image, int *histogram);
	int static getThreshold(Image &image);
	void static toGrayScale(Image &image);
	int static Binarize(Image &image, int threshold = -1);
	void static Incline(Image &image, double R,double G,double B);
	void static Sobel(Image &image);
	void static Laplacian(Image &image);
	void static Gaussian(Image &image, float sigma);
	PointList static Harris(Image &image, float sigma = 1, float threshold = 0.1);
	void static AntiNoise(Image &image, unsigned int level = 1);
	void static Contrast(Image &image);
	void static YIQRange(Image &image, int Ylow, int Yhigh,int Ilow = 0, int Ihigh = 0, int Qlow = 0, int Qhigh = 0);
	void static HSVRange(Image &image, int Hlow, int Hhigh,int Slow = 0, int Shigh = 100, int Vlow = 0, int Vhigh = 100);
	void static HRange(Image &image, int low, int high, bool in_range = true); 
	void static SRange(Image &image, int low, int high, bool in_range = true); 
	void static VRange(Image &image, int low, int high, bool in_range = true); 
	void static GammaCorrection(Image &image, double gamma);
	void static Brightness(Image &image, int threshold=128);
	Point static TemplateMatching(Image &image, Image temp, int threshold = -1);
	Point static TemplateMatching(Image &image, int x, int y, int width, int height, Image temp, int threshold = -1);
	int static Matching(Image &image, int x,int y, Image temp);
	PointList static getNodeList(Image &image,int connected_num);
	int static getConnectedNum(Image &image,int x,int y);
	int static Labeling(Image &image);
	int static SamplingObject(Image &image, int label);
	int static SamplingLargeObject(Image &image);
	double static getBulr(Image &image);
	
};
