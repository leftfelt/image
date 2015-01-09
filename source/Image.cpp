#include <Image.hpp>
#include <Math.hpp>
#include <Pixel.hpp>
#include <functional>

void Image::for_each(Image &image, std::function<void(int x,int y)> func){
	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			func(i,j);
		}
	}
}

void Image::around_each(Image &image, Point point, int around, std::function<void(int x, int y)> func){
	for(int j = point.Y()-around ; j <= point.Y()+around ; j++){
		for(int i = point.X()-around ; i <= point.X()+around ; i++){
			func(i,j);
		}
	}
}



Image::Image(){
	this->init();
}

Image::Image(int width, int height){
	this->init();
	this->Create(width,height);
}

Image::Image(const Image &image){
	mimage = image.mimage;
	this->clearColor = image.clearColor;
	this->width =image.width;
	this->height = image.height;
}

Image::~Image(){
}

void Image::init(){
	clearColor = Pixel();

	height = 0;
	width = 0;

	mimage.clear();
	std::vector<Pixel>(mimage).swap(mimage);
}

void Image::Create(int width, int height){
	this->init();

	this->width = width;
	this->height = height;

	
	mimage.resize((ulong)this->width * this->height);
}

void Image::Paste(int x, int y, Image &image, int mix_rate){
	Image::for_each(image,[&](int i, int j){
		if(image.Get(i,j) == Pixel() || i%(int)(1.0/(mix_rate/100.0))) return;
		Pixel pixel = image.Get(i,j);
		this->Put((x+i),(y+j),pixel);
	});
}


Image Image::Cut(int x, int y, int width, int height){
	Image ret(width,height);
	Image::for_each(ret,[&](int i, int j){
		ret.Put(i,j,this->Get((x+i),(y+j)));
	});

	return ret;
}

void Image::Delete(){
	this->init();
}


Pixel Image::Clear(){
	return this->clearColor;
}


void Image::Clear(Pixel color){
	this->clearColor = color;
}


int Image::Width(){
	return this->width;
}

int Image::Height(){
	return this->height;
}


void Image::Width(int width){
	Size(width,this->Height());
}

void Image::Height(int height){
	Size(this->Width(),height);
}

void Image::Size(int width ,int height ){
	Image image;

	double wRate = (double)this->Width() / width;
	double hRate = (double)this->Height() / height;

	width = abs(width);
	height = abs(height);

	image.Create(width,height);

	Image::for_each(image,[&](int i,int j){
		image.Put(i,j,this->Get((int)(i*wRate),(int)(j*hRate)));
	});
	image.Clear(this->clearColor);
	*this = image;
}

Image Image::Rotate(int angle){
	Image image;
	int rx,ry;
	double rad = angle * (3.14159265358979323846 / 180.0);
	
	image.Create(Width(),Height());
	image = this->clearColor;
	image.Clear(this->clearColor);

	for(int j = -this->Height()/2 ; j <= this->Height()/2 ; j++){
		for(int i = -this->Width()/2 ; i <= this->Width()/2 ; i++){
			rx = (int)(this->Width()/2 + i * cos(rad) - j * sin(rad));
			ry = (int)(this->Height()/2 + i * sin(rad) + j * cos(rad));
			if(rx < 0 || this->Width() <= rx || ry < 0 || this->Height() <= ry)continue;
			image.Put(i+this->Width()/2,j+this->Height()/2,Get(rx,ry));
		}
	}
	return image;
}


Image Image::operator=(Pixel pixel){
	
	
	std::fill(this->mimage.begin(), this->mimage.end(), pixel);
	
	return  *this;
}


Image Image::operator+(Image image){
	Image temp;

	temp = *this;
	temp += image;
	
	return temp;
}

Image Image::operator-(Image image){
	Image temp;

	temp = *this;
	temp -= image;
	
	return temp;
}
Image Image::operator*(Image image){
	Image temp;

	temp = *this;
	temp *= image;
	
	return temp;
}
Image Image::operator/(Image image){
	Image temp;

	temp = *this;
	temp /= image;
	
	return temp;
}

Image Image::operator+(Pixel pixel){
	Image temp;

	temp = *this;
	temp += pixel;
	
	return temp;
}

Image Image::operator-(Pixel pixel){
	Image temp;

	temp = *this;
	temp -= pixel;
	
	return temp;
}
Image Image::operator*(Pixel pixel){
	Image temp;

	temp = *this;
	temp *= pixel;
	
	return temp;
}
Image Image::operator/(Pixel pixel){
	Image temp;

	temp = *this;
	temp /= pixel;
	
	return temp;
}

Image Image::operator+=(Image image){	
	int width,height;

	
	width = (this->Width() > image.Width()) ? image.Width() : this->Width();
	height = (this->Height() > image.Height()) ? image.Height() : this->Height();
	
	for(int j = 0 ; j < height ; j++){
		for(int i = 0 ; i < width ; i++){
			this->mimage.at((ulong)i+j*width) = this->mimage.at((ulong)i+j*width) + image.mimage.at((ulong)i+j*width);
		}
	}

	return *this;
}


Image Image::operator-=(Image image){	
	int width,height;

	
	width = (this->Width() > image.Width()) ? image.Width() : this->Width();
	height = (this->Height() > image.Height()) ? image.Height() : this->Height();

	for(int j = 0 ; j < height ; j++){
		for(int i = 0 ; i < width ; i++){
			this->mimage.at((ulong)i+j*width) = this->mimage.at((ulong)i+j*width) - image.mimage.at((ulong)i+j*width);
		}
	}

	return *this;
}

Image Image::operator*=(Image image){	
	int width,height;

	
	width = (this->Width() > image.Width()) ? image.Width() : this->Width();
	height = (this->Height() > image.Height()) ? image.Height() : this->Height();

	for(int j = 0 ; j < height ; j++){
		for(int i = 0 ; i < width ; i++){
			this->mimage.at((ulong)i+j*width) = this->mimage.at((ulong)i+j*width) * image.mimage.at((ulong)i+j*width);
		}
	}

	return *this;
}


Image Image::operator/=(Image image){	
	int width,height;

	
	width = (this->Width() > image.Width()) ? image.Width() : this->Width();
	height = (this->Height() > image.Height()) ? image.Height() : this->Height();

	for(int j = 0 ; j < height ; j++){
		for(int i = 0 ; i < width ; i++){
			this->mimage.at((ulong)i+j*width) = this->mimage.at((ulong)i+j*width) / image.mimage.at((ulong)i+j*width);
		}
	}

	return *this;
}

Image Image::operator+=(Pixel pixel){	
	Image::for_each(*this,[&](int i, int j){
		this->mimage.at((ulong)i+j*this->Width()) = this->mimage.at((ulong)i+j*this->Width()) + pixel;
	});
	return *this;
}

Image Image::operator-=(Pixel pixel){	
	Image::for_each(*this,[&](int i, int j){
		this->mimage.at((ulong)i+j*this->Width()) = this->mimage.at((ulong)i+j*this->Width()) - pixel;
	});
	return *this;
}

Image Image::operator*=(Pixel pixel){	
	Image::for_each(*this,[&](int i, int j){
		this->mimage.at((ulong)i+j*this->Width()) = this->mimage.at((ulong)i+j*this->Width()) * pixel;
	});
	return *this;
}

Image Image::operator/=(Pixel pixel){	
	Image::for_each(*this,[&](int i, int j){
		this->mimage.at((ulong)i+j*this->Width()) = this->mimage.at((ulong)i+j*this->Width()) / pixel;
	});
	return *this;
}