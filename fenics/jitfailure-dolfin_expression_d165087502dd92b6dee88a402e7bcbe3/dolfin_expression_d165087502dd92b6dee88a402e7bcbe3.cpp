
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
  class dolfin_expression_d165087502dd92b6dee88a402e7bcbe3 : public Expression
  {
     public:
       double viscosity;
double wp;
double ww;
double r0;
double angle1;
double angle2;
double k;
double d;


       dolfin_expression_d165087502dd92b6dee88a402e7bcbe3()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = 6*viscosity*((wp*(x[0]*r0+d)+ww*r0*x[0])*sin(angle1)-((ww+wp)*r0*x[1])*sin(angle2))/k-duiliu;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "viscosity") { viscosity = _value; return; }          if (name == "wp") { wp = _value; return; }          if (name == "ww") { ww = _value; return; }          if (name == "r0") { r0 = _value; return; }          if (name == "angle1") { angle1 = _value; return; }          if (name == "angle2") { angle2 = _value; return; }          if (name == "k") { k = _value; return; }          if (name == "d") { d = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "viscosity") return viscosity;          if (name == "wp") return wp;          if (name == "ww") return ww;          if (name == "r0") return r0;          if (name == "angle1") return angle1;          if (name == "angle2") return angle2;          if (name == "k") return k;          if (name == "d") return d;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_d165087502dd92b6dee88a402e7bcbe3()
{
  return new dolfin::dolfin_expression_d165087502dd92b6dee88a402e7bcbe3;
}

