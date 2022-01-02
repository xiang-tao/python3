
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
  class dolfin_expression_577f4be50a397044c0ccdf6a7f558b69 : public Expression
  {
     public:
       std::shared_ptr<dolfin::GenericFunction> generic_function_ppp;
std::shared_ptr<dolfin::GenericFunction> generic_function_CF;
double dH;
std::shared_ptr<dolfin::GenericFunction> generic_function_p;
double px_p;


       dolfin_expression_577f4be50a397044c0ccdf6a7f558b69()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double p;
            generic_function_p->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&p), x);
          double CF;
            generic_function_CF->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&CF), x);
          double ppp;
            generic_function_ppp->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&ppp), x);
          values[0] = r0/p0*(ppp*CF*dH+3*p*p*CF*px_p);

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "dH") { dH = _value; return; }          if (name == "px_p") { px_p = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "dH") return dH;          if (name == "px_p") return px_p;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "ppp") { generic_function_ppp = _value; return; }          if (name == "CF") { generic_function_CF = _value; return; }          if (name == "p") { generic_function_p = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "ppp") return generic_function_ppp;          if (name == "CF") return generic_function_CF;          if (name == "p") return generic_function_p;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_577f4be50a397044c0ccdf6a7f558b69()
{
  return new dolfin::dolfin_expression_577f4be50a397044c0ccdf6a7f558b69;
}

