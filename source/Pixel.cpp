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

unsigned char Pixel::Red(){
	return this->red;
}
unsigned char Pixel::Green(){
	return this->green;
}
unsigned char Pixel::Blue(){
	return this->blue;
}
unsigned char Pixel::Lightness(){
	return (unsigned char)(( this->Red() + this->Green() + this->Blue() ) / 3);
}
int Pixel::Label(){
	return this->label;
}

void Pixel::Red(unsigned char red){
	this->empty = false;
	this->red = red;
}
void Pixel::Green(unsigned char green){
	this->empty = false;
	this->green = green;
}
void Pixel::Blue(unsigned char blue){
	this->empty = false;
	this->blue = blue;
}
void Pixel::Lightness(unsigned char lightness){
	this->empty = false;
	this->Red(lightness);
	this->Green(lightness);
	this->Blue(lightness);
}
void Pixel::Label(int label){
	this->empty = false;
	this->label = label;
}


unsigned char Pixel::Y(){
	return (unsigned char)( this->Red() * 0.299 + this->Green() * 0.587 + this->Blue() * 0.114);
}
char Pixel::I(){
	return (char)( this->Red() * 0.596 - this->Green() * 0.274 - this->Blue() * 0.322);
}
char Pixel::Q(){
	return (char)( this->Red() * 0.211 - this->Green() * 0.522 - this->Blue() * 0.311);
}


int Pixel::H(){
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

unsigned char Pixel::S(){
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

unsigned char Pixel::V(){
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


Pixel Pixel::operator+(Pixel pixel){
	this->Red( this->Red() + pixel.Red() );
	this->Green( this->Green() + pixel.Green() );
	this->Blue( this->Blue() + pixel.Blue() );
	return *this;
}
Pixel Pixel::operator-(Pixel pixel){
	this->Red( this->Red() - pixel.Red() );
	this->Green( this->Green() - pixel.Green() );
	this->Blue( this->Blue() - pixel.Blue() );
	return *this;
}
Pixel Pixel::operator*(Pixel pixel){
	this->Red( this->Red() * pixel.Red() );
	this->Green( this->Green() * pixel.Green() );
	this->Blue( this->Blue() * pixel.Blue() );
	return *this;
}
Pixel Pixel::operator/(Pixel pixel){
	if(pixel.Red() == 0) this->Red(0);
	else this->Red( this->Red() / pixel.Red() );

	if(pixel.Green() == 0) this->Green(0);
	else this->Green( this->Green() / pixel.Green() );

	if(pixel.Blue() == 0) this->Blue(0);
	else this->Blue( this->Blue() / pixel.Blue() );
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


Pixel Pixel::operator+(unsigned char value){
	this->Red(this->Red() + value);
	this->Green(this->Green() + value);
	this->Blue(this->Blue() + value);
	return *this;
}
Pixel Pixel::operator-(unsigned char value){
	this->Red(this->Red() - value);
	this->Green(this->Green() - value);
	this->Blue(this->Blue() - value);
	return *this;
}
Pixel Pixel::operator*(unsigned char value){
	this->Red(this->Red() * value);
	this->Green(this->Green() * value);
	this->Blue(this->Blue() * value);
	return *this;
}
Pixel Pixel::operator/(unsigned char value){
	this->Red(this->Red() / value);
	this->Green(this->Green() / value);
	this->Blue(this->Blue() / value);
	return *this;
}

Pixel Pixel::operator+=(unsigned char value){
	*this = *this + value;
	return *this;
}
Pixel Pixel::operator-=(unsigned char value){
	*this = *this - value;
	return *this;
}
Pixel Pixel::operator*=(unsigned char value){
	*this = *this * value;
	return *this;
}
Pixel Pixel::operator/=(unsigned char value){
	*this = *this / value;
	return *this;
}



bool Pixel::operator==(Pixel pixel) const {
	return (this->empty == pixel.empty && (
				this->red	 == pixel.Red()
				&& this->green == pixel.Green()
				&& this->blue	 == pixel.Blue()
			));
}

bool Pixel::operator!=(Pixel pixel){
	return !(*this == pixel);
}