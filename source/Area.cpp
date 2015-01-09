#include <Area.hpp>


Area::Area(int x=0, int y=0, int width=0, int height=0){
	mx			= x;
	my			= y;
	mwidth		= width;
	mheight	= height;
}

Area::Area(){
	mx			= 0;
	my			= 0;
	mwidth		= 0;
	mheight	= 0;
}


void Area::Pos(int x, int y){
	mx = x;
	my = y;
}


void Area::Size(int width, int height){
	mwidth = width;
	mheight = height;
}


int Area::X(){
	return mx;
}


int Area::Y(){
	return my;
}


int Area::Width(){
	return mwidth;
}


int Area::Height(){
	return mheight;
}


bool Area::HitCheck(int x,int y){
	if( mx <= x && x <= mx+mwidth && my <= y && y <= my+mheight){
		return true;
	}
	return false;
}

bool Area::HitCheck(Area area){
	if( area.HitCheck(mx,my) ||
		area.HitCheck(my+mwidth,my) ||
		area.HitCheck(my,my+mheight) ||
		area.HitCheck(my+mwidth,my+mheight) 
		){
		return true;
	}
	return false;
}
