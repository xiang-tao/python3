
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_3e54374ae35206c4fcbed3ff8634596c : public Expression
  {
     public:
       double H0;
double xx;
double angle1;
std::shared_ptr<dolfin::GenericFunction> generic_function_Hhat;
std::shared_ptr<dolfin::GenericFunction> generic_function_Hhat_dx;
std::shared_ptr<dolfin::GenericFunction> generic_function_Hhat_dxx;
std::shared_ptr<dolfin::GenericFunction> generic_function_p;


       dolfin_expression_3e54374ae35206c4fcbed3ff8634596c()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double p;
            generic_function_p->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&p), x);
          double Hhat_dxx;
            generic_function_Hhat_dxx->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&Hhat_dxx), x);
          double Hhat_dx;
            generic_function_Hhat_dx->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&Hhat_dx), x);
          double Hhat;
            generic_function_Hhat->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&Hhat), x);
          values[0] = 3*H0*H0/p0*p*p*Hhat*(-xx*sin(angle1))*Hhat_dx+H0*H0/p0*p*p*p*(Hhat_dx*Hhat_dx+Hhat*Hhat_dxx);

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "H0") { H0 = _value; return; }          if (name == "xx") { xx = _value; return; }          if (name == "angle1") { angle1 = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "H0") return H0;          if (name == "xx") return xx;          if (name == "angle1") return angle1;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "Hhat") { generic_function_Hhat = _value; return; }          if (name == "Hhat_dx") { generic_function_Hhat_dx = _value; return; }          if (name == "Hhat_dxx") { generic_function_Hhat_dxx = _value; return; }          if (name == "p") { generic_function_p = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "Hhat") return generic_function_Hhat;          if (name == "Hhat_dx") return generic_function_Hhat_dx;          if (name == "Hhat_dxx") return generic_function_Hhat_dxx;          if (name == "p") return generic_function_p;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_3e54374ae35206c4fcbed3ff8634596c()
{
  return new dolfin::dolfin_expression_3e54374ae35206c4fcbed3ff8634596c;
}

