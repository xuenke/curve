#ifndef __Fit_h__
#define __Fit_h__
 
#include <vector>
#include <math.h>
 
template<int degrees>
class CFit
{
public:
  CFit(std::vector<double>& data_points_x, std::vector<double>& data_points_y):data_points_x_(data_points_x), 
                                                                               data_points_y_(data_points_y), 
                                                                               ssr_(0.0), sse_(0.0), rmse_(0.0) {
    b_fit_yet_ = false;
  }
 
  ~CFit() {}
protected: 
  // gaussian elimination
  template<int real_degree>
  static void GaussianElimination(double (&matrix)[real_degree+1][real_degree+2], 
                                  double (&coe_array)[real_degree+1]) {
    int i, j, k;
    for (i = 0; i < real_degree; i++ ) {          // loop to perform the gauss elimination
      for (k = i + 1; k < (real_degree+1); k++) {
        double t = matrix[k][i]/matrix[i][i];
        for (j = 0; j <= (real_degree+1); j++)
            matrix[k][j] -= t*matrix[i][j];     // make the elements below the pivot elements equal to zero or elimnate the variables
      }
    }
    // back-substitution
    for (i = real_degree; i >= 0; i--) {        //x is an array whose values correspond to the values of x,y,z..
      coe_array[i] = matrix[i][real_degree+1];  //make the variable to be calculated equal to the rhs of the last equation
      for (j = 0; j < (real_degree+1); j++)
        if (j != i)                               //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
            coe_array[i] -= matrix[i][j]*coe_array[j];
      coe_array[i] = coe_array[i] / matrix[i][i];        //now finally divide the rhs by the coefficient of the variable to be calculated
    }
  }
 
  ///
  /// \brief calculate y value according to approximation function
  /// \return y value respecte to x
  ///
  template<typename T>
  double GetY(const T x) const {
    double ans(0);
    for (int i = 0; i < (degrees+1); ++i) {
      ans += coefficients_[i] * pow((double)x, (int)i);
    }
    return ans;
  }
 
  ///
  /// \brief calculate mean value
  /// \return mean value
  ///
  template <typename T>
  static T Mean(const std::vector<T>& v) {
    return Mean(&v[0], v.size());
  }

  template <typename T>
  static T Mean(const T* v, int length) {
    T total(0);
    for (int i = 0; i < length; ++i) {
      total += v[i];
    }
    return (total / length);
  }
 
  template<typename T>
  void CalcError(const T* x, const T* y, int length,
    double& r_ssr, double& r_sse, double& r_rmse) {
    T mean_y = Mean<T>(y, length);
    T yi(0);
      
    for (int i = 0; i < length; ++i) {
      yi = GetY(x[i]);
      r_ssr += ((yi-mean_y) * (yi-mean_y)); // explained Sum of Squares
      r_sse += ((yi-y[i]) * (yi-y[i])); // residual sum of squares
    }
    r_rmse = sqrt(r_sse / (double(length)));
  }
 
 
  /**
  * @brief  implement polynomial approximation
  * @author 
  * @param  [in]  int points_num, count of data points
            [in]  const std::vector<double>& data_points_x, x value of data points
            [in]  const std::vector<double>& data_points_y, y value of data points
  * @param  [out]   double (&coefficientArr)[degrees+1], coefficients of the approximated polynomial, from low to high
  * @return   none
  * @note   none
  */
  static void PolynomialFit(int points_num, const std::vector<double>& data_points_x,
                            const std::vector<double>& data_points_y, double (&coefficientArr)[degrees+1]) {
    int i = 0, j = 0, k = 0;
    double x[2*degrees+1] = {0};          //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    for (i = 0; i < 2 * degrees + 1; i++) {
      for (j = 0; j < points_num; j++)
        x[i] += pow(data_points_x[j], i);          //consecutive positions of the array will store points_num,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    }
    double y[degrees+1] = {0};          //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    for (i = 0; i < degrees + 1; i++) {   
      for (j = 0; j < points_num; j++)
        y[i] += pow(data_points_x[j], i) * data_points_y[j];   //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    }

    double equation_matrix[degrees+1][degrees+2] = {0};       //equation_matrix is the Normal matrix(augmented) that will store the equations
    for (i = 0; i <= degrees; i++) {
      for (j = 0; j <= degrees; j++) {
        equation_matrix[i][j] = x[i+j];             //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
      }
      equation_matrix[i][degrees+1] = y[i];          //load the values of y as the last column of equation_matrix(Normal Matrix but augmented)
    }

    GaussianElimination<degrees>(equation_matrix, coefficientArr);
  }

public:
  // polynomial approximation interface
  void PolyFit() {
    if (data_points_x_.size() == data_points_y_.size()) {
      PolynomialFit(static_cast<int>(data_points_x_.size()), data_points_x_, data_points_y_, coefficients_);
      b_fit_yet_ = true;
      CalcError(&data_points_x_[0], &data_points_y_[0], static_cast<int>(data_points_x_.size()), ssr_, sse_, rmse_);
    }
  }

  //- calculate y value accoding to x
  double UnaryPolynomialCalc(double dx) {
    double dy = 0.0;
    for (int i = 0; i <= degrees; ++i) {
      dy += pow(dx, (double)i) * coefficients_[i];
    }
    return b_fit_yet_ ? dy : 0.0;
  }

  // get a series of points from polynomail curve according to given resolution and length
  void Frenet2Cartisian(double start_x, double path_len, 
                  std::vector<double>& p_x, std::vector<double>& p_y) {
    p_x.clear();
    p_y.clear();
    double resolution = 0.15;
    double dx = 0.0;
    for (int i = 1; i*resolution < path_len; i++) {
      start_x += dx;
      p_x.push_back(start_x);
      p_y.push_back(UnaryPolynomialCalc(start_x));
      double derivative = DerivativeCalc(start_x);
      double alpha = atan2(derivative, 1.0);
      dx = resolution * cos(alpha);
    }
  }

  // calculate derivative of the polynomial at x
  double DerivativeCalc(double x) {
    double d = 0.0;
    for (int i = 1; i <= degrees; i++) {
      d += pow(x, double(i-1)) * coefficients_[i] * i;
    }
    return d;
  }

  ///
  /// \brief calculate residual sum of squares
  /// \return residual sum of squares
  ///
  double GetSSE() {return sse_;}
  ///
  /// \brief explained sum of squares
  /// \return explained sum of squares
  ///
  double GetSSR() {return ssr_;}
  ///
  /// \brief mean squared error
  /// \return mean squared error
  ///
  double GetRMSE() {return rmse_;}
  ///
  /// \brief goodness-of-fit
  /// \return goodness-of-fit
  ///
  double GetGoodnessOfFit() {return 1 - (sse_ / (ssr_ + sse_));}
 
  ///
  /// \brief get coefficient according to order
  /// \return coefficient of a certatin order
  ///
  double GetFactor (int i)
  {
    return (i <= degrees) ? coefficients_[i] : 0.0;
  }
   
private:
  double coefficients_[degrees+1];
  const std::vector<double>& data_points_x_;
  const std::vector<double>& data_points_y_;
  bool b_fit_yet_;// whether the fitting process have done yet
 
  double ssr_;           // explained sum of square
  double sse_;           // residual sum of square
  double rmse_;          // root mean square error
};
 
 
#endif // __Fit_h__
