#pragma once

#include <math.h>

#define AreaList std::vector<Area>
/*
矩形又は円形のエリアを作成する
エリア同士の衝突判定を行う
*/
class Area{
private:
	int mx;			
	int my;			
	int mheight;	
	int mwidth;		
public:
	Area();
	Area(int x, int y, int width, int height);		
	void Pos(int x, int y);					
	void Size(int width, int height);		
	int X();							
	int Y();							
	int Width();						
	int Height();						
	bool HitCheck(int x,int y);				
	bool HitCheck(Area area);				
};
