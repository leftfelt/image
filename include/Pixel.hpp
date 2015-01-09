#ifndef PIXELCLASS_INCLUDE
#define PIXELCLASS_INCLUDE

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

	Pixel operator+(int value);
	Pixel operator-(int value);
	Pixel operator*(int value);
	Pixel operator/(int value);
	Pixel operator+=(int value);
	Pixel operator-=(int value);
	Pixel operator*=(int value);
	Pixel operator/=(int value);

	bool operator==(Pixel pixel);
	bool operator!=(Pixel pixel);
};


inline unsigned char Pixel::Red(){
	return this->red;
}
inline unsigned char Pixel::Green(){
	return this->green;
}
inline unsigned char Pixel::Blue(){
	return this->blue;
}
inline unsigned char Pixel::Lightness(){
	return ( this->Red() + this->Green() + this->Blue() ) / 3;
}
inline int Pixel::Label(){
	return this->label;
}


inline void Pixel::Red(unsigned char red){
	this->empty = false;
	this->red = red;
}
inline void Pixel::Green(unsigned char green){
	this->empty = false;
	this->green = green;
}
inline void Pixel::Blue(unsigned char blue){
	this->empty = false;
	this->blue = blue;
}
inline void Pixel::Lightness(unsigned char lightness){
	this->empty = false;
	this->Red(lightness);
	this->Green(lightness);
	this->Blue(lightness);
}
inline void Pixel::Label(int label){
	this->empty = false;
	this->label = label;
}


inline unsigned char Pixel::Y(){
	return (unsigned char)( this->Red() * 0.299 + this->Green() * 0.587 + this->Blue() * 0.114);
}
inline char Pixel::I(){
	return (char)( this->Red() * 0.596 - this->Green() * 0.274 - this->Blue() * 0.322);
}
inline char Pixel::Q(){
	return (char)( this->Red() * 0.211 - this->Green() * 0.522 - this->Blue() * 0.311);
}


inline int Pixel::H(){
	unsigned char R,G,B,RGB[3];
	unsigned char max=0,min=255;
	double h = 0;
	
	R = this->Red();
	G = this->Green();
	B = this->Blue();

	RGB[0] = R;
	RGB[1] = G;
	RGB[2] = B;
	
	int i;
	for(i = 0 ; i < 3 ; i++){
		if(max < RGB[i])max = RGB[i];
		if(min > RGB[i])min = RGB[i];
	}

	if( max == R )	h = (60 * ( (float)( G-B )/( max-min )+0 ));
	if( max == G )	h = (60 * ( (float)( B-R )/( max-min )+2 ));
	if( max == B )	h = (60 * ( (float)( R-G )/( max-min )+4 ));
	if( max == 0 ) h = 0;

	if( h < 0) h += 360;

	return (int)h;
}

inline unsigned char Pixel::S(){
	unsigned char RGB[3];
	double max=0.0,min=1.0;
	int i;
	unsigned char s = 0;

	RGB[0] = this->Red();
	RGB[1] = this->Green();
	RGB[2] = this->Blue();

	for(i = 0 ; i < 3 ; i++){
		if(max < RGB[i]/255.0)max = RGB[i]/255.0;
		if(min > RGB[i]/255.0)min = RGB[i]/255.0;
	}

	s =  (unsigned char)((max - min)/max * 100.0);

	if( max == 0) s = 0;

	return s;
}

inline unsigned char Pixel::V(){
	unsigned char RGB[3];
	double max=0.0;
	int i;

	RGB[0] = this->Red();
	RGB[1] = this->Green();
	RGB[2] = this->Blue();

	for(i = 0 ; i < 3 ; i++){
		if(max < RGB[i]/255.0)max = RGB[i]/255.0;
	}
	return (unsigned char)(max * 100);
}

#endif