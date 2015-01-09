#include <ImageUtil.hpp>
#include <string.h>

float ImageUtil::laplacian_filter[9] = {-1,-1,-1, -1, 8,-1, -1,-1,-1};

float ImageUtil::Gaussian(int x, int y, float sigma){
	return 1.0/(2.0*Math::Pi*sigma) * exp( -(pow((float)x, 2)+pow((float)y, 2))/(2.0*pow((float)sigma, 2)) );
}


float ImageUtil::derivativesGaussian(int x, int y, float sigma){
	return (float)(-x)/(2.0*Math::Pi*pow((float)sigma, 4)) * exp( -(pow((float)x, 2)+pow((float)y, 2))/(2.0*pow((float)sigma, 2)) );
}

double ImageUtil::getBulr(Image &image){
	double num = 0;
	Pixel pixel;
	Image::for_each(image,[&](int i, int j){
		pixel = image.Get(i,j);
		num += (pixel.Red() + pixel.Blue() + pixel.Green());
	});
	num /= (double)(image.Width() * image.Height());
	return (double)(( 1.0 / num ) * 100);
}

void ImageUtil::Shrink(Image &image){
	int cnt;
	Image output(image.Width(),image.Height());

	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			cnt = 0;
			
			if(image.Get(i,j).Lightness() == 0){
				output.Put(i,j,image.Get(i,j));
				continue;
			}
		
			
			for(int k = -1 ; k <= 1  ; k++){
				for(int l = -1 ; l <= 1 ; l++){
					if(i+l<0 || j+k<0 || i+l>=image.Width() || j+k>=image.Height())continue;
					if(k != 1 || l != 1){
						if(image.Get(i+l, j+k).Lightness() == 0){
							cnt++;
						}
					}
				}
			}

			if(cnt != 0){
				output.Put(i,j,Pixel(0));
			}else{
				output.Put(i,j,image.Get(i,j));
			}
		}
	}

	image = output;
	output.Delete();
}
void ImageUtil::Expand(Image &image){
	int cnt;
	Image output(image.Width(),image.Height());

	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			cnt = 0;
			
			if(image.Get(i,j).Lightness() == 255){
				output.Put(i,j,image.Get(i,j));
				continue;
			}
		
			
			for(int k = -1 ; k <= 1  ; k++){
				for(int l = -1 ; l <= 1 ; l++){
					if(i+l<0 || j+k<0 || i+l>=image.Width() || j+k>=image.Height())continue;
					if(k != 1 || l != 1){
						if(image.Get(i+l, j+k).Lightness() == 255){
							cnt++;
						}
					}
				}
			}

			if(cnt != 0){
				output.Put(i,j,Pixel(255));
			}else{
				output.Put(i,j,image.Get(i,j));
			}
		}
	}

	image = output;
	output.Delete();
}

void ImageUtil::Opening(Image &image,int level){
	
	for(int i = 0 ; i < level ; i++){
		Shrink(image);
	}

	
	for(int i = 0 ; i < level ; i++){
		Expand(image);
	}
}


void ImageUtil::Closing(Image &image,int level){
	
	for(int i = 0 ; i < level ; i++){
		Expand(image);
	}

	
	for(int i = 0 ; i < level ; i++){
		Shrink(image);
	}
}


void ImageUtil::TopHat(Image &image,int level){
	Image temp;
	temp = image;
	Opening(temp,level);
	image = image - temp;

}

void ImageUtil::BlackHat(Image &image,int level){
	Image temp;
	temp = image;
	Closing(temp,level);
	image = image - temp;
}

int ImageUtil::Thinning(Image &image){

	Image temp;
	
	int mask[][9] = {
		{ 0, 0,-1,
		  0, 1, 1,
		 -1, 1,-1},
		{ 0, 0, 0,
		 -1, 1,-1,
		  1, 1,-1},
		{-1, 0, 0,
		  1, 1, 0,
		 -1, 1,-1},
		{ 1,-1, 0,
		  1, 1, 0,
		 -1,-1, 0},
		{-1, 1,-1,
		  1, 1, 0,
		 -1, 0, 0},
		{-1, 1, 1,
		 -1, 1,-1,
		  0, 0, 0},
		{-1, 1,-1,
		  0, 1, 1,
		  0, 0,-1},
		{ 0,-1,-1,
		  0, 1, 1,
		  0,-1, 1},
		{-2},
	};

	
	int count = 0;
	temp = image;
	bool continue_flag = true;
	
	int mask_num;
	while(continue_flag){
		continue_flag = false;
		mask_num = 0;
		while(mask[mask_num][0] != -2){
			for(int j=1 ; j<image.Height()-1 ; j++){
				for(int i=1 ; i<image.Width()-1 ; i++){
					if(image.Get(i,j) == Pixel(0)) continue;
					int k;
					for(k=0 ; k<9 ; k++){
						
						if(mask[mask_num][k] != -1){
							int pixel = mask[mask_num][k] * 255;
							
							if(image.Get(i+(k%3)-1,j+(k/3)-1) != Pixel(pixel))break;
						}
					}
					
					if(k == 9){
						
						temp.Put(i,j,Pixel(0));
						count++;
						continue_flag = true;
					}
				}
			}
			image = temp;
			mask_num++;
		}
	}
	return count;
}


void ImageUtil::UnSharpMasking(Image &image,int level){
	Image temp;
	temp = image;
	LowResolution(temp,level);
	image = image + (image - temp)*2;
}


void ImageUtil::LowResolution(Image &image,int level){
	
	long R=0,G=0,B=0;
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			R=0;
			G=0;
			B=0;
			for(int k=-level ; k<=level ; k++){
				for(int l=-level ; l<=level ; l++){
					Pixel pixel = image.Get(i+l,j+k);
					R += pixel.Red();
					G += pixel.Green();
					B += pixel.Blue();
				}
			}
			image.Put(i,j,Pixel(
				(unsigned char)(R/pow((level*2+1)*1.0,2)),
				(unsigned char)(G/pow((level*2+1)*1.0,2)),
				(unsigned char)(B/pow((level*2+1)*1.0,2))
			));
		}
	}
}


void ImageUtil::colorRegionSplit(Image &image,int threshold){
	int x,y;
	int label_num = 0;
	int unattended_num = 0;
	std::map<unsigned long,unsigned long> label_disp; 
	std::map<unsigned long,Pixel> tab; 
	std::map<unsigned long,unsigned long> label_tab_num;
	std::map<unsigned long,int> unlabel_x;
	std::map<unsigned long,int> unlabel_y;

	for(y = 0 ; y < image.Height() ; y++){
		for(x = 0 ; x < image.Width() ; x++){
			if(label_disp[x+y*image.Width()] > 0) continue;
			label_num++;
			label_disp[x+y*image.Width()] = label_num;
			tab[label_num] = image.Get(x,y);
			unattended_num = 1;
			unlabel_x[unattended_num]=x; 
			unlabel_y[unattended_num]=y;
			while(unattended_num > 0){
				int mx = unlabel_x[unattended_num];
				int my = unlabel_y[unattended_num];
				unattended_num--;
				for(int i=0 ; i<9 ; i++){
					int cx = (i%3)-1;
					int cy = (i/3)-1;
					if(0<=mx+cx && mx+cx<image.Width() &&
							0<=my+cy && my+cy<image.Height()){
						if(label_disp[mx+cx+(my+cy)*image.Width()] == 0){
							Pixel pixel = image.Get(mx+cx,my+cy);
							int diff =  abs(tab[label_num].Red()-pixel.Red()) +
										abs(tab[label_num].Green()-pixel.Green()) +
										abs(tab[label_num].Blue()-pixel.Blue());
							if(diff > threshold) continue;
							label_disp[mx+cx+(my+cy)*image.Width()] = label_num;
							label_tab_num[label_num]++;
							unattended_num++;
							unlabel_x[unattended_num] = mx+cx;
							unlabel_y[unattended_num] = my+cy;
						}
					}
				}
			}
		}
	}
	int seg_num = label_num; 
	for(y = 0 ; y < image.Height() ; y++){
		for(x = 0 ; x < image.Width() ; x++){
			label_num = label_disp[x+y*image.Width()];
			Pixel pixel = tab[label_num];
			image.Put(x,y,pixel);
		}
	}
}

void ImageUtil::getHistogram(Image &image, int *histogram){
	int i;
	for(i = 0; i < 256 ; i++){
		*(histogram+i)=0;
	}
	for(int j = 0 ; j<image.Height() ; j++){
		for(int i = 0; i < image.Width() ; i++){
			(*(histogram+image.Get(i,j).Lightness()))++;
		}
	}
}


int ImageUtil::getThreshold(Image &image){
	int histogram[256];
	int i,j;
	int a,b;
	int sum = 0;
	int n1 = 0,n2 = 0;
	double u0,u1,u2;
	double O12 = 0,O22 = 0;
	double OW2;
	double OB2;
	int t = 0;
	int threshold = 0;
	double proportion = -1;

	getHistogram(image,histogram);
	
	
	for(i = 0; i <= 255 ; i++){
		if((*(histogram+i)) != 0){
			a = i;
			break;
		}
	}
	
	for(i = 255; i >= 0 ; i--){
		if((*(histogram+i)) != 0){
			b = i;
			break;
		}
	}
	
	for(i = 0 ; i <= 255 ; i++){
		sum += i*(*(histogram+i));
	}

	
	u0 = sum / (image.Width()*image.Height());

	for(t = a ; t <= b ; t++){
		
		sum = 0;
		n1 = 0;
		for(i = a ; i <= t ; i++){
			sum += i*(*(histogram+i));
			n1 += *(histogram+i);
		}

		u1 = sum / n1; 
	
		
		O12 = 0;
		for(i = a ; i <= t ; i++){
			for(j = 0 ; j < *(histogram+i) ; j++){
				O12 += pow((i - u0),2);
			}
		}

		O12 = O12 / n1;
		sum = 0;
		n2 = 0;
		for(i = t ; i <= b ; i++){
			sum += i*(*(histogram+i));
			n2 += *(histogram+i);
		}

		u2 = sum / n2; 

		O22 = 0;
		
		for(i = t ; i <= b ; i++){
			for(j = 0 ; j < *(histogram+i) ; j++){
				O22 += pow((i - u0),2);
			}
		}

		O22 = O22 / n2;

		OW2 = ( n1 * O12 + n2 * O22 ) / (n1 + n2);
		OB2 = ( n1*pow((u1 - u0),2) + n2*pow((u2-u0),2) ) / (n1+n2);

		if(proportion < (OB2/OW2)){
			proportion = OB2/OW2;
			threshold = t;
		}
	}

	return threshold;
}

int ImageUtil::Binarize(Image &image, int threshold){
	int ret = 0;
	
	if(threshold == -1){
		ret = getThreshold(image);
		threshold = ret;
	}
	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			if(image.Get(i,j).Lightness() >= threshold){
				image.Put(i,j,Pixel(255));
			}else{
				image.Put(i,j,Pixel(0));
			}
		}
	}

	return ret;

}


void ImageUtil::toGrayScale(Image &image){
	Image::for_each(image,[&](int x, int y){
		image.Put(x, y, image.Get(x, y).Lightness());
	});
}


void ImageUtil::Incline(Image &image, double R,double G,double B){
	
	
	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			image.Put(
				i,j,
				Pixel(
					(int)(R * image.Get(i,j).Red()),
					(int)(G * image.Get(i,j).Green()),
					(int)(B * image.Get(i,j).Blue())
				)
			);
		}
	}
}

void ImageUtil::filtering(Image &image, float *filter, int filter_size[2], int direction[2]){
	Image temp = image;
	int i, j, height = image.Height(), width = image.Width();
	for(j = 0 ; j < height ; j++){
		for(i = 0 ; i < width ; i++){
			temp.Put(i,j,fabs(ImageUtil::convolution(image, filter, i, j, filter_size, direction)));
		}
	}
	image = temp;
}

float ImageUtil::convolution(Image &image, float *filter, int x, int y, int filter_size[2], int direction[2]){
	float pixel = 0;
	int i, ix, iy, iwidth = filter_size[0]*direction[0]*direction[1];
	
	for(i = 0 ; i < filter_size[0]*direction[0] ; i++){
		ix = i*direction[0];
		iy = i*direction[1];
		pixel += (float)image.Get(x+ix,y+iy).Lightness() * filter[i+i*iwidth];
	}
	
	for(i = 0 ; i < filter_size[1]*direction[1] ; i++){
		ix = i*direction[1];
		iy = i*direction[0];
		pixel += (float)image.Get(x+ix,y+iy).Lightness() * filter[i+i*iwidth];
	}

	return pixel;
}




void ImageUtil::Sobel(Image &image){
	int width, height, graylevel;
	int pos = 0;
	int num[9];
	int dx,dy,gray;
	Image outimage;

	outimage = image;
	width = image.Width();
	height = image.Height();
	graylevel = 255;

	float filter_y[9] = {1, 2, 1, 0, 0, 0, -1, -2, -1};
	float filter_x[9] = {1, 0, -1, 2, 0, -2, 1, 0, -1};
	int filter_size[2] = {3, 3};
	int direction_x[2] = {1,1};
	int direction_y[2] = {1,1};

	
	
	
	
	for(int j = 0 ; j < height ; j++){
		for(int i = 0 ; i < width ; i++){
			
			dx = ImageUtil::convolution(image, filter_x, i, j, filter_size,direction_x);
			dy = ImageUtil::convolution(image, filter_y, i, j, filter_size,direction_y);

			gray = abs(dx)+abs(dy);

			if(gray > graylevel) gray = graylevel;
			outimage.Put(i,j,Pixel(gray));
		}
	}
	image = outimage;
}


void ImageUtil::Laplacian(Image &image){
	int filter_size[2] = {3, 3};
	int direction[2] = {1,1};
	ImageUtil::filtering(image, ImageUtil::laplacian_filter, filter_size, direction);
}


void ImageUtil::Gaussian(Image &image, float sigma){
	float gaussian_filter[9];
	int ix = 0, iy = 0;
	int filter_size[2] = {3, 3};
	int direction[2] = {1,1};

	for(int i = 0 ; i < 9 ; i++){
		ix = i%3-1;
		iy = i/3-1;
		gaussian_filter[i] = ImageUtil::Gaussian(ix, iy, sigma);
	}

	ImageUtil::filtering(image, gaussian_filter,filter_size, direction);
}

PointList ImageUtil::Harris(Image &image, float sigma, float threshold){
	PointList point_list;
	float gaussian_x[3], gaussian[9];
	int ix = 0, iy = 0;
	int filter_size[2] = {3, 3};
	int direction[2] = {1,1};
	int direction_x[2] = {1,0};
	int direction_y[2] = {0,1};

	Image imx = image, imy = image;
	
	for(int i = 0 ; i < 3 ; i++){
		gaussian_x[i] = ImageUtil::derivativesGaussian(i-1, 0, sigma);
	}

	for(int i = 0 ; i < 9 ; i++){
		ix = i%3-1;
		iy = i/3-1;
		gaussian[i] = ImageUtil::Gaussian(ix, iy, sigma);
	}

	ImageUtil::filtering(imx, gaussian_x,filter_size, direction_x);
	ImageUtil::filtering(imy, gaussian_x,filter_size, direction_y);

	
	Image Wxx = imx * imx, Wxy = imx * imy, Wyy = imy * imy;
	ImageUtil::filtering(Wxx, gaussian,filter_size,direction);
	ImageUtil::filtering(Wxy, gaussian,filter_size,direction);
	ImageUtil::filtering(Wyy, gaussian,filter_size,direction);

	
	Image Wdet = Wxx*Wyy - Wxy*2, Wtr = Wxx + Wyy;

	Image::for_each(Wdet,[&](int x, int y){
		float pixel = (((float)Wdet.Get(x,y).Lightness() / (float)Wtr.Get(x,y).Lightness()));
		if(threshold <= pixel){
			image.Put(x,y, 255);
			point_list.push_back(Point(x,y));
		}
	});

	
	return point_list;
}


void ImageUtil::AntiNoise(Image &image,unsigned int level){
	Opening(image,level);
}


void ImageUtil::Contrast(Image &image){
	int histogram[256];
	int a = 0,b = 255;

	getHistogram(image,histogram);

	
	for(int i = 0; i <= 255 ; i++){
		if(histogram[i] != 0){
			a = i;
			break;
		}
	}

	for(int i = 255; i >= 0 ; i--){
		if(histogram[i] != 0){
			b = i;
			break;
		}
	}

	
	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			if( image.Get(i,j).Lightness() >= 0 && image.Get(i,j).Lightness() < a){
				image.Put(i,j,Pixel(0));
			}else if(image.Get(i,j).Lightness() >=a && image.Get(i,j).Lightness() <= b){
				image.Put(i,j,Pixel( (int)(255 * ((image.Get(i,j).Lightness() - a) / (double)(b - a))) ));
			}else if(image.Get(i,j).Lightness() > b && image.Get(i,j).Lightness() <= 255){
				image.Put(i,j,Pixel(255));
			}
		}
	}
}


void ImageUtil::YIQRange(Image &image, int Ylow, int Yhigh,int Ilow, int Ihigh, int Qlow, int Qhigh){
	int yiq_I,yiq_Q,yiq_Y;

	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			yiq_Y = image.Get(i,j).Y();
			yiq_I = image.Get(i,j).I();
			yiq_Q = image.Get(i,j).Q();
			if( (Ylow <= yiq_Y && yiq_Y <= Yhigh) &&
				((Ilow <= yiq_I && yiq_I <= Ihigh) || (Ilow == Ihigh)) && 
				((Qlow <= yiq_Q && yiq_Q <= Qhigh) || (Qlow == Qhigh)) 
				){ 
					image.Put(i,j,Pixel(255));
			}else{
					image.Put(i,j,Pixel(0));
			}
		}
	}
}


void ImageUtil::HSVRange(Image &image, int Hlow, int Hhigh,int Slow, int Shigh, int Vlow, int Vhigh){
	ImageUtil::HRange(image,Hlow,Hhigh);
	ImageUtil::SRange(image,Slow,Shigh);
	ImageUtil::VRange(image,Vlow,Vhigh);
}


void ImageUtil::HRange(Image &image, int low, int high, bool in_range){
	Pixel pixel;
	low		= Math::limit(low	,0,360);
	high	= Math::limit(high	,0,360);
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			pixel = image.Get(i,j);
			if( !((low <= pixel.H() && pixel.H() <= high) == in_range) ) {
				image.Put(i,j,Pixel(0));
			}
		}
	}
}

void ImageUtil::SRange(Image &image, int low, int high, bool in_range){
	Pixel pixel;
	low		= Math::limit(low	,0,100);
	high	= Math::limit(high	,0,100);
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			pixel = image.Get(i,j);
			if( !((low <= pixel.S() && pixel.S() <= high) == in_range) ) {
				image.Put(i,j,Pixel(0));
			}
		}
	}
}

void ImageUtil::VRange(Image &image, int low, int high, bool in_range){
	Pixel pixel;
	low		= Math::limit(low	,0,100);
	high	= Math::limit(high	,0,100);
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			pixel = image.Get(i,j);
			if( !((low <= pixel.V() && pixel.V() <= high) == in_range) ) {
				image.Put(i,j,Pixel(0));
			}
		}
	}
}

void ImageUtil::GammaCorrection(Image &image, double gamma){
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			Pixel pixel = image.Get(i,j);
			image.Put(
				i,j,Pixel(
					(unsigned char)(255*pow((pixel.Red()/255.0),1/gamma)),
					(unsigned char)(255*pow((pixel.Green()/255.0),1/gamma)),
					(unsigned char)(255*pow((pixel.Blue()/255.0),1/gamma))
				));
		}
	}
}


void ImageUtil::Brightness(Image &image, int threshold){
	int average = 0;
	int count = 0;
	for(int j = 0 ; j < image.Height() ; j++){
		for(int i = 0 ; i < image.Width() ; i++){
			Pixel pixel = image.Get(i,j);
			if(pixel.Lightness() != 255){
				average += pixel.Lightness();
				count++;
			}
		}
	}
	average /= count;
	image += (threshold - average);
}




Point ImageUtil::TemplateMatching(Image &image, Image temp, int threshold){
	int i,j,k,l;
	Point pt(0,0);
	int max = 255 * 3 *temp.Width() * temp.Height();
	int min = max;
	int valid_num = 0;
	int count = 0;

	
	if(image.Width() < temp.Width() || 
		image.Height() < temp.Height())return pt;

	for(j = -temp.Height() ; j < image.Height() ; j++){
		for(i = -temp.Width() ; i < image.Width() ; i++){
			count = 0;
			valid_num = 0;
			for(l = 0 ; l < temp.Height(); l++){
				for(k = 0 ; k < temp.Width() ; k++){
					
					if( temp.Get(k,l) == temp.Clear()) {
						count += 255*3;
						continue;
					}
					if((i+k) >= image.Width() || (j+l) >= image.Height() || (i+k) < 0 || (j+l) < 0){
						count += 255*3;
						continue;
					}
					valid_num++;
					count += abs( image.Get((i+k),(j+l)).Blue()		- temp.Get(k,l).Blue());
					count += abs( image.Get((i+k),(j+l)).Green()	- temp.Get(k,l).Green());
					count += abs( image.Get((i+k),(j+l)).Red()		- temp.Get(k,l).Red());
				}
			}
			max = 255*3*valid_num;
			
			if(min > count && (threshold <= (max-count)/(double)max * 100 || threshold == -1) ){
				min = count;
				pt.Pos(i,j);
			}
		}
	}

	return pt;
}


Point ImageUtil::TemplateMatching(Image &image, int x, int y, int width, int height, Image temp, int threshold){
	int i,j;
	Point pt(x,y);
	int max = 0;
	int count = 0;

	
	if(image.Width() < temp.Width() || 
		image.Height() < temp.Height())return pt;

	for(j = y-height/2 ; j <= y+height/2 ; j++){
		for(i = x-width/2 ; i <= x+width/2 ; i++){
			count = 0;
			if(i >= image.Width() || j >= image.Height() || i < 0 || j < 0)continue;
			count = Matching(image,i-temp.Width()/2,j-temp.Height()/2,temp);
			
			if(max < count && (threshold <= count || threshold == -1) ){ 
				max = count;
				pt.Pos(i,j);
			}
		}
	}

	return pt;
}


int ImageUtil::Matching(Image &image, int x,int y, Image temp){
	int i,j;
	int max = 255 * 3 *temp.Width() * temp.Height();
	int count = 0;
	int valid_num = 0;

	
	if(image.Width() < temp.Width() || 
		image.Height() < temp.Height())return 0;

	for(j = 0 ; j < temp.Height(); j++){
		for(i = 0 ; i < temp.Width() ; i++){
			
			if( temp.Get(i,j) == temp.Clear()){
				
				continue;
			}
			if((i+x) >= image.Width() || (j+y) >= image.Height() || (i+x) < 0 || (j+y) < 0){
				
				count += 255*3;
				continue;
			}
			valid_num++;
			count += abs( image.Get((i+x),(j+y)).Blue()  - temp.Get(i,j).Blue());
			count += abs( image.Get((i+x),(j+y)).Green() - temp.Get(i,j).Green());
			count += abs( image.Get((i+x),(j+y)).Red()   - temp.Get(i,j).Red());
		}
	}
	max = 255*3*valid_num;
	return (int)((max-count)/(double)max * 100);
}


PointList ImageUtil::getNodeList(Image &image,int connected_num){
	PointList point_list;
	Point point;
	for(int j = 0 ; j < image.Height(); j++){
		for(int i = 0 ; i < image.Width() ; i++){
			if(image.Get(i,j) == Pixel(0)) continue;
			if(ImageUtil::getConnectedNum(image,i,j) == connected_num){
				point.Pos(i,j);
				point_list.push_back(point);
			}
		}
	}
	return point_list;
}


int ImageUtil::getConnectedNum(Image &image,int x,int y){
	int count = 0;
	for(int j = -1 ; j <= 1; j++){
		for(int i = -1 ; i <= 1 ; i++){
			if(i==0&&j==0)continue;
			if(image.Get(x+i,y+j) == Pixel(0)) continue;
			count++;
		}
	}
	return count;
}



int ImageUtil::Labeling(Image &image){
	int w,h,j,k,i,size;
	int label_count=1,flag_count;
	int min_data[8],min;
	int *Lookup_table,*sub_table,*label_list;
	int max=0,label=0,Lookup_flag=0;
	
	Lookup_table = new int[image.Width() * image.Height()];
	sub_table = new int[image.Width() * image.Height()];
	label_list = new int[image.Width() * image.Height()];

	
	size = image.Width() * image.Height();
	for(i = 0 ;i< size;i++){
		*(label_list+i) = 0;
		*(sub_table+i) = 0;
	}
	for(i=0;i<size;i++){
		Lookup_table[i] = i;
	}

	
	
	for(k=0;k<2;k++){
		for(h = 1; h < image.Height()-1;h++){
			for(w=1;w<image.Width()-1;w++){
				if(image.Get(w,h).Lightness() == 255){
					flag_count = 0;
					for(int i=0 ; i<9 ; i++){
						if( (i%3-1) == 0 && (i/3-1) == 0 ) continue;
						if( *(label_list+(w+i%3-1)+(h+i/3-1)*image.Width()) == 0 ) flag_count++;
					}
					
					if(flag_count == 8){
						
					    
						
						*(label_list+w+h*image.Width()) = label_count;
						label_count++;
					}else{
						
					    
						
						for(int i=0,count=0 ; i<9 ; i++){
							if( (i%3-1) == 0 && (i/3-1) == 0 ) continue;
							min_data[count] = *(label_list+(w+i%3-1)+(h+i/3-1)*image.Width());
							count++;
						}
	
						min = min_calc(min_data);
						
					    
						
	
						for(j=0;j<8;j++){
							if(min_data[j] != min && min_data[j] != 0){
								if(Lookup_table[min_data[j]] >= min){
									Lookup_table[min_data[j]] = min;
								}
							}
						}
											
						*(label_list+w+h*image.Width()) = min;
					}
				}
			}
		}
		
		for(h=0;h<label_count;h++){
			if(h != Lookup_table[h]){
				Lookup_table[h] = Lookup_update(h,Lookup_table);
			}
		}

			
		for(h = 0 ; h < image.Height() ; h++){
			for(w = 0 ; w < image.Width() ; w++){
				if(*(label_list+w+h*image.Width())  != 0){
					*(label_list+w+h*image.Width())  = Lookup_table[*(label_list+w+h*image.Width())];
				}
			}
		}
	}

	
	i=0;
	sub_table = (int *)memset(sub_table,-1,sizeof(int) * image.Width() * image.Height());
	for(h=0;h<label_count;h++){
		if(sub_table[Lookup_table[h]] == -1){
			sub_table[Lookup_table[h]] = i;
			i++;
		}
	}
	label = i;
	
	for(h=0;h<label_count;h++){
		Lookup_table[h] = sub_table[Lookup_table[h]];
	}

	Pixel pixel;
	
	for(h = 0 ; h < image.Height() ; h++){
		for(w = 0 ; w < image.Width() ; w++){
			if(*(label_list+w+h*image.Width()) != 0){
				pixel = image.Get(w,h);
				pixel.Label(Lookup_table[*(label_list+w+h*image.Width())]);
				image.Put(w,h,pixel);
			}
		}
	}
	
	free(Lookup_table);
	free(sub_table);

	return label;
}


int ImageUtil::SamplingObject(Image &image, int label){
	Pixel pixel;
	int count=0;
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			pixel = image.Get(i,j);
			if((pixel.Label() == label)){
				pixel.Lightness(255);
				image.Put(i,j,pixel);
				count++;
			}else{
				pixel.Lightness(0);
				image.Put(i,j,pixel);
			}
		}
	}

	return count;
}


int ImageUtil::SamplingLargeObject(Image &image){
	std::map<int,int> area_map;
	std::map<int,int>::iterator area;
	
	for(int j=0 ; j < image.Height() ; j++){
		for(int i=0 ; i < image.Width() ; i++){
			int id = image.Get(i,j).Label();
			if(id != 0){
				area_map[id]++;
			}
		}
	}
	
	int area_id, num,max=0;
	for(area = area_map.begin(); area != area_map.end() ; ++area){
 		num = area->second;
		if(max < num){
			max = num;
			area_id = area->first;
		}
	}
	ImageUtil::SamplingObject(image,area_id);
	return area_id;
}



int ImageUtil::min_calc(int *p)
{
	int i,min=0xFFFF;

	for(i=0;i<8;i++){
		if( *(p+i)!=0 ){
				if(min > *(p+i)){
				min = *(p+i);
			}
		}
	}

	return min;
}

int ImageUtil::Lookup_update(int i,int * label)
{
	int a;

	if(i == label[i]){
		return label[i];
	}else{
		a =  Lookup_update(label[i],label);
	}

	return a;
}