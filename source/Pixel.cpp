#include <Pixel.hpp>
#include <Math.hpp>


Pixel::Pixel(){
	this->initialize();
}
Pixel::Pixel(unsigned char lightness){
	this->initialize();
	this->Lightness(lightness);
}
Pixel::Pixel(unsigned char red,unsigned char green,unsigned char blue){
	this->initialize();
	this->Red(red);
	this->Green(green);
	this->Blue(blue);
}


Pixel::~Pixel(){
}

void Pixel::setRGB(unsigned char red,unsigned char green,unsigned char blue){
	this->Red(red);
	this->Green(green);
	this->Blue(blue);
}

void Pixel::initialize(){
	this->red = 0;
	this->green = 0;
	this->blue = 0;
	this->label = 0;
	this->empty = true;
}





Pixel Pixel::operator+(Pixel pixel){
	this->Red( Math::limit(this->Red() + pixel.Red(),0,255) );
	this->Green( Math::limit(this->Green() + pixel.Green(),0,255) );
	this->Blue( Math::limit(this->Blue() + pixel.Blue(),0,255) );
	return *this;
}
Pixel Pixel::operator-(Pixel pixel){
	this->Red( Math::limit(this->Red() - pixel.Red(),0,255) );
	this->Green( Math::limit(this->Green() - pixel.Green(),0,255) );
	this->Blue( Math::limit(this->Blue() - pixel.Blue(),0,255) );
	return *this;
}
Pixel Pixel::operator*(Pixel pixel){
	this->Red( Math::limit(this->Red() * pixel.Red(),0,255) );
	this->Green( Math::limit(this->Green() * pixel.Green(),0,255) );
	this->Blue( Math::limit(this->Blue() * pixel.Blue(),0,255) );
	return *this;
}
Pixel Pixel::operator/(Pixel pixel){
	if(pixel.Red() == 0) this->Red(0);
	else this->Red( Math::limit(this->Red() / pixel.Red(),0,255) );

	if(pixel.Green() == 0) this->Green(0);
	else this->Green( Math::limit(this->Green() / pixel.Green(),0,255) );

	if(pixel.Blue() == 0) this->Blue(0);
	else this->Blue( Math::limit(this->Blue() / pixel.Blue(),0,255) );
	return *this;
}

Pixel Pixel::operator+=(Pixel pixel){
	*this = *this + pixel;
	return *this;
}
Pixel Pixel::operator-=(Pixel pixel){
	*this = *this - pixel;
	return *this;
}
Pixel Pixel::operator*=(Pixel pixel){
	*this = *this * pixel;
	return *this;
}
Pixel Pixel::operator/=(Pixel pixel){
	*this = *this / pixel;
	return *this;
}


Pixel Pixel::operator+(int value){
	this->Red(  Math::limit(this->Red() + value,0,255) );
	this->Green(  Math::limit(this->Green() + value,0,255) );
	this->Blue(  Math::limit(this->Blue() + value,0,255) );
	return *this;
}
Pixel Pixel::operator-(int value){
	this->Red(  Math::limit(this->Red() - value,0,255) );
	this->Green(  Math::limit(this->Green() - value,0,255) );
	this->Blue(  Math::limit(this->Blue() - value,0,255) );
	return *this;
}
Pixel Pixel::operator*(int value){
	this->Red(  Math::limit(this->Red() * value,0,255) );
	this->Green(  Math::limit(this->Green() * value,0,255) );
	this->Blue(  Math::limit(this->Blue() * value,0,255) );
	return *this;
}
Pixel Pixel::operator/(int value){
	this->Red(  Math::limit(this->Red() / value,0,255) );
	this->Green(  Math::limit(this->Green() / value,0,255) );
	this->Blue(  Math::limit(this->Blue() / value,0,255) );
	return *this;
}

Pixel Pixel::operator+=(int value){
	*this = *this + value;
	return *this;
}
Pixel Pixel::operator-=(int value){
	*this = *this - value;
	return *this;
}
Pixel Pixel::operator*=(int value){
	*this = *this * value;
	return *this;
}
Pixel Pixel::operator/=(int value){
	*this = *this / value;
	return *this;
}



bool Pixel::operator==(Pixel pixel){
	return (this->empty == pixel.empty && (
				this->Red()	 == pixel.Red()
				&& this->Green() == pixel.Green()
				&& this->Blue()	 == pixel.Blue()
			));
}

bool Pixel::operator!=(Pixel pixel){
	return !(*this == pixel);
}