#ifndef POINTCLASS_INCLUDE
#define POINTCLASS_INCLUDE

#define PointList std::vector<Point>
/*
*/
class Point{
private:
	int mx;			
	int my;			
public:
	Point();
	Point(int x, int y);		
	void Pos(int x, int y);				
	int X();							
	int Y();							

	Point operator+(Point &point);
	Point operator-(Point &point);
	Point operator*(Point &point);
	Point operator/(Point &point);
	
	Point operator=(int value);
	Point operator+(int value);
	Point operator-(int value);
	Point operator*(int value);
	Point operator/(int value);
	
	Point operator+=(int value);
	Point operator-=(int value);
	Point operator*=(int value);
	Point operator/=(int value);
};

#endif