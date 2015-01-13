#include <ImageIO.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

bool ImageIO::Save(Image &image, std::string filename){
	std::ofstream out;

	out.open(filename);
	if(!out) {
		return false;
	}
	out << "P3\n";
	out << image.Width() << " " << image.Height() << "\n";
	out << "255\n";
	Image::for_each(image,[&](int x, int y){
		Pixel pixel = image.Get(x,y);
		out << (int)pixel.Red() << " " << (int)pixel.Green() << " " << (int)pixel.Blue() << std::endl;
	});
	out.close();
	return true;
}

bool ImageIO::Load(Image &image, std::string filename){
	return ImageIO::LoadBitmap(image, filename);
}


bool ImageIO::LoadPPM(Image &image, std::string filename){
	//*
	std::ifstream in;
	std::string data;
	int width, height;
	unsigned char r,g,b;
	char *ptr = 0;

	in.open(filename);
	if(!in) {
		return false;
	}

	while(getline(in, data)) if(data[0] != '#') break;
	
	
	while(in >> data){if(data[0] != '#'){width = (int)std::strtol(data.c_str(),&ptr,10);break;}}
	while(in >> data){if(data[0] != '#'){height = (int)std::strtol(data.c_str(),&ptr,10);break;}}

	while(in >> data) if(data[0] != '#') break;

	image.Create(width, height);
	
	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			while(in >> data){if(data[0] != '#'){r = (unsigned char)std::strtol(data.c_str(),&ptr,10);break;}}
			while(in >> data){if(data[0] != '#'){g = (unsigned char)std::strtol(data.c_str(),&ptr,10);break;}}
			while(in >> data){if(data[0] != '#'){b = (unsigned char)std::strtol(data.c_str(),&ptr,10);break;}}
			image.Put(x, y, Pixel(r, g, b));
		}
	}

	in.close();
	//*/
}

bool ImageIO::LoadBitmap(Image &image, std::string filename){
	FILE *fp;

	int width, height;
	short bitCount ;
	int data_ofset;
	unsigned char R,G,B;
	
	fp = fopen( filename.c_str(),"rb");
	if(fp == NULL ){
		return false;
	}

	
	image.Delete();

	
	fseek(fp,sizeof(unsigned char)*10,SEEK_CUR);
	fread( &data_ofset, sizeof (unsigned int),1, fp );
	fseek(fp,sizeof(unsigned char)*4,SEEK_CUR);
	fread( &width, sizeof(long),1,fp );
	fread( &height, sizeof(long),1,fp);
	fseek(fp,sizeof(unsigned char)*2,SEEK_CUR);
	fread( &bitCount, sizeof(unsigned short),1,fp);

	fseek(fp,data_ofset,SEEK_SET);

	image.Create(width, height);

	
	for(int j=0 ; j<image.Height() ; j++){
		for(int i=0 ; i<image.Width() ; i++){
			if(bitCount == 8){
				fread( &G	, sizeof(unsigned char), 1, fp );
				image.Put(i,j,Pixel(G));
			}else{
				fread( &B	, sizeof(unsigned char), 1, fp );
				fread( &G	, sizeof(unsigned char), 1, fp );
				fread( &R	, sizeof(unsigned char), 1, fp );

				image.Put(i,j,Pixel(R,G,B));

				if(bitCount == 32)fread( &G, sizeof (unsigned char), 1, fp );
			}
		}
		
		if(bitCount == 24){
			for(int k = 0 ; k < image.Width()%4 ; k++){
				fread( &G, sizeof (unsigned char), 1, fp );
			}
		}
	}

	
	image.Height(-image.Height());

	fclose(fp);

	return true;
}