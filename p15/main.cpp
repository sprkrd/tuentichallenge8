#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

struct Lcg {
  uint128_t u, i, o; // multiplier, increment, modulus
  uint128_t operator()(uint128_t x) const {
    return (x*u+i)%o;
  }
};

uint128_t firstRepetition(uint128_t s, const Lcg& f) {
  uint128_t t,h,mu,lambda;
  t = f(s);
  h = f(f(s));
  while (t != h) {
    t = f(t);
    h = f(f(h));
  }
  mu = 0;
  t = s;
  while (t != h) {
    t = f(t);
    h = f(h);
    ++mu;
  }
  lambda = 1;
  h = f(h);
  while (t != h) {
    h = f(h);
    ++lambda;
  }
  return mu+lambda+1;
}

int main() {
  int number_of_cases;
  cin >> number_of_cases;
  for (int case_i = 1; case_i <= number_of_cases; ++case_i) {
    uint128_t s,u,i,o;
    cin >> s >> u >> i >> o;
    auto first_repetition = firstRepetition(s,Lcg{u,i,o});
    cout << "Case #" << case_i << ": " << first_repetition << endl;
  }
}
