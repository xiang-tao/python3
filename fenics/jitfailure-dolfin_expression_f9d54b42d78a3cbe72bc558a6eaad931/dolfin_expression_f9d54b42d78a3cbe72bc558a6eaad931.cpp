
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
  class dolfin_expression_f9d54b42d78a3cbe72bc558a6eaad931 : public Expression
  {
     public:
       std::shared_ptr<dolfin::GenericFunction> generic_function_ppp;
std::shared_ptr<dolfin::GenericFunction> generic_function_p;
double r0;
std::shared_ptr<dolfin::GenericFunction> generic_function_g;
double py_p;
double p0;
double ww;


       dolfin_expression_f9d54b42d78a3cbe72bc558a6eaad931()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double g;
            generic_function_g->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&g), x);
          double p;
            generic_function_p->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&p), x);
          double ppp;
            generic_function_ppp->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&ppp), x);
          values[0] = r0*r0/p0*(ww*ww*ppp+3*g*p*p/r0*py_p;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "r0") { r0 = _value; return; }          if (name == "py_p") { py_p = _value; return; }          if (name == "p0") { p0 = _value; return; }          if (name == "ww") { ww = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "r0") return r0;          if (name == "py_p") return py_p;          if (name == "p0") return p0;          if (name == "ww") return ww;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "ppp") { generic_function_ppp = _value; return; }          if (name == "p") { generic_function_p = _value; return; }          if (name == "g") { generic_function_g = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "ppp") return generic_function_ppp;          if (name == "p") return generic_function_p;          if (name == "g") return generic_function_g;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_f9d54b42d78a3cbe72bc558a6eaad931()
{
  return new dolfin::dolfin_expression_f9d54b42d78a3cbe72bc558a6eaad931;
}

