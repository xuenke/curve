#include "polynomial_fit.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

struct Vector3d
{
  double x;
  double y;
  double z;
};


bool writeAsc( const std::string& filename, const std::vector<Vector3d>& points )
{
	std::ofstream fout( filename.c_str());
	if( fout.fail() )
		return false;

	for( int i = 0; i!= points.size(); ++i )
	{
		fout << points[i].x << "," << points[i].y << "," << points[i].z << std::endl;
	}
	fout.close();

	return true;

}

void CreateRoughPath(std::vector<double>& p_x, std::vector<double>& p_y,
        std::vector<double>& p_theta) {
    double radius = 5.0;
    double max_peak = 0.1;
    double max_bias = 0.1;
    int num = 50;
    double step = M_PI / num / 2.0;
    p_x.clear();
    p_y.clear();
    p_theta.clear();
    srand(time(NULL));
    for (int i = num - 1; i >= 0; i--) {
        double rand_y = rand() % 100 / 100.0 * max_peak;
        double rand_theta = rand() % 100 / 100.0 * max_bias + M_PI;
        p_x.push_back(radius * cos(i * step));
        p_y.push_back(radius * sin(i * step)+ rand_y);
        p_theta.push_back(M_PI / 2 + i * step + rand_theta);
    }
}

void Cartisian2Frenet() {

}

int main() {
  std::vector<Vector3d> points;
  // double y[] = {7,16,6,18,6,6,10,8};
  // double x[] = {-109.71,-101.81,-103.83,-99.89,-90,-112.17,-93.5,-96.13};
  // double x[] = {0.0, 5.0, 7.0, 8.0, 10.0, 15.0, 18.0};
  // double y[] = {3.0, -2.0, -5.0, 0.0, 6.0, 12.0, 8.0};
  // double x[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
  // double y[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
  // std::vector<double> xArr(std::begin(x),std::end(x));
  // std::vector<double> yArr(std::begin(y),std::end(y));

  std::vector<double> xArr;
  std::vector<double> yArr;
  std::vector<double> thetaArr;
  CreateRoughPath(xArr, yArr, thetaArr);

  // std::vector<double> xArr, yArr, dArr;
  // for (int i = 0; i < 10; i++) {
  //   xArr.push_back(double(i));
  //   yArr.push_back(double(1.0 + 2.0 * i + 3.0 * i*i + 4.0 * i*i*i));
  //   dArr.push_back(double(2.0 + 6.0 * i + 12.0 * i*i));
  //   Vector3d point;
  //   point.x = xArr.back();
  //   point.y = yArr.back();
  //   point.z = dArr.back();
  //   points.push_back(point);
  // }
  // writeAsc("origin.csv", points);

  typedef CFit<4> LineFit;
  LineFit objPolyfit(xArr,yArr);
  objPolyfit.PolyFit();

  std::cout << "fitting done....." << std::endl;

  std::vector<double> p_x, p_y;
  objPolyfit.Frenet2Cartisian(0.0, 10.0, p_x, p_y);

  points.clear();
  for (int i = 0; i < xArr.size(); i++) {
    Vector3d point;
    point.x = xArr[i];
    point.y = yArr[i];
    point.z = objPolyfit.DerivativeCalc(xArr[i]);
    points.push_back(point);
  }
  writeAsc("fit.csv", points);

  points.clear();
  for (int i = 0; i < p_x.size(); i++) {
    Vector3d point;
    point.x = p_x[i];
    point.y = p_y[i];
    point.z = 0.0;
    points.push_back(point);
  }
  writeAsc("fit2.csv", points);

  return 0;
}