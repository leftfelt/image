#pragma once

#include <Math.h>
/*
* PixelClass
*/

class Pixel{
private:
	unsigned char red;
	unsigned char green;
	unsigned char blue;
	int label;
	bool empty;
	void initialize();
public:
	Pixel();
	Pixel(unsigned char lightness);
	Pixel(unsigned char red,unsigned char green,unsigned char blue);
	~Pixel();
	void setRGB(unsigned char red,unsigned char green,unsigned char blue);
	unsigned char Red();
	unsigned char Green();
	unsigned char Blue();
	unsigned char Lightness();
	int Label();
	void Red(unsigned char red);
	void Green(unsigned char green);
	void Blue(unsigned char blue);
	void Lightness(unsigned char lightness);
	void Label(int label);
	
	unsigned char Y();
	char I();
	char Q();
	
	int H();
	unsigned char S();
	unsigned char V();
	
	Pixel operator+(Pixel pixel);
	Pixel operator-(Pixel pixel);
	Pixel operator*(Pixel pixel);
	Pixel operator/(Pixel pixel);
	Pixel operator+=(Pixel pixel);
	Pixel operator-=(Pixel pixel);
	Pixel operator*=(Pixel pixel);
	Pixel operator/=(Pixel pixel);

	Pixel operator+(unsigned char value);
	Pixel operator-(unsigned char value);
	Pixel operator*(unsigned char value);
	Pixel operator/(unsigned char value);
	Pixel operator+=(unsigned char value);
	Pixel operator-=(unsigned char value);
	Pixel operator*=(unsigned char value);
	Pixel operator/=(unsigned char value);

	bool operator==(Pixel pixel);
	bool operator!=(Pixel pixel);
};