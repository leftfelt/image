#include <Math.hpp>

double Math::Pi = 3.14159265358979323846264338327950288 ; 

double Math::Angle(double radian){
	return (double)radian / (Pi / 180);
}
double Math::Radian(double angle){
	return (double)angle * ( Pi / 180);
}

double Math::distance(Point point1, Point point2){
	return sqrt( pow(point2.X() - point1.X(),2.0) + pow(point2.Y() - point1.Y(),2.0) );
}


double Math::generateUniform(std::mt19937 rng){
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(rng);
}

double Math::generateRandomNumber(double min, double max, std::mt19937 rng){
	if ((double)0 <= min || max < (double)0) return min + (Math::generateUniform(rng)) * (max - min);
	return Math::generateUniform(rng) * ((rand()%2) ? max : min);
}
