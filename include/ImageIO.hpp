#pragma once

#include <Image.hpp>
#include <string>

class ImageIO {
private:
	bool static LoadBitmap(Image &image, std::string filename);
public:
	bool static Load(Image &image, std::string filename);
	bool static LoadPPM(Image &image, std::string filename);
	bool static Save(Image &image, std::string filename);
};
