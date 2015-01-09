#ifndef BIIMAGECLASS_INCLUDE
#define BIIMAGECLASS_INCLUDE

#include <Pixel.hpp>
#include <Point.hpp>
#include <functional>
#include <vector>

class Image {
protected:
	std::vector<Pixel> mimage;
private:
	Pixel clearColor;
	int width;
	int height;
	void init();
public:
	void static for_each(Image &image, std::function<void(int x,int y)> func);
	void static around_each(Image &image, Point point, int around, std::function<void(int x, int y)> func);
	Image();
	Image(int width, int height);
	Image(const Image &image);
	~Image();
	void Create(int width, int height);
	void Delete();
	void Paste(int x, int y, Image &image, int mix_rate = 100);
	Image Cut(int x, int y, int width, int height);
	Image Rotate(int angle);
	
	Pixel Clear();
	void Clear(Pixel color);
	int Width();
	void Width(int width);
	int Height();
	void Height(int height);
	void Size(int width,int height);
	void Put(int x, int y, Pixel pixel);
	Pixel Get(int x,int y);
	
	Image operator=(Pixel pixel);
	Image operator+(Image image);
	Image operator-(Image image);
	Image operator*(Image image);
	Image operator/(Image image);
	Image operator+(Pixel pixel);
	Image operator-(Pixel pixel);
	Image operator*(Pixel pixel);
	Image operator/(Pixel pixel);
	Image operator+=(Image image);
	Image operator-=(Image image);
	Image operator*=(Image image);
	Image operator/=(Image image);
	Image operator+=(Pixel pixel);
	Image operator-=(Pixel pixel);
	Image operator*=(Pixel pixel);
	Image operator/=(Pixel pixel);
};

inline void Image::Put(int x, int y, Pixel pixel){
	x = ((x < 0)? this->Width()+x : x)%this->Width();
	y = ((y < 0)? this->Height()+y : y)%this->Height();
	mimage.at(x+y*this->Width()) = pixel;
}

inline Pixel Image::Get(int x,int y){
	x = ((x < 0)? this->Width()+x : x)%this->Width();
	y = ((y < 0)? this->Height()+y : y)%this->Height();
	if( this->clearColor == mimage.at(x+y*this->Width()) ) return Pixel();
	return mimage.at(x+y*this->Width());
}

#endif