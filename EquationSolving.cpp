#include <iostream>
#include <vector>
#include <iomanip>
namespace EqSol {
	const double eps = 1e-14;
	double f(double x)
	{
		return x * x - 2;
	}
	double binsolution(const double a, const double b)
	{
		double l = a, r = b;
		if (abs(f(l)) < eps)
			return l;
		if (abs(f(r)) < eps)
			return r;
		while (r - l > eps)
		{
			double m = (l + r) / 2;
			if (signbit(f(a)) != signbit(f(m)))
				r = m;
			else l = m;
		}
		return r;
	}
}
int main()
{
	double first = EqSol::binsolution(-2, 0);
	double second = EqSol::binsolution(0, 2);
	std::cout << std::setprecision(13) <<first <<"      Diverge: "<<abs(first+sqrt(2))<<'\n';
	std::cout << std::setprecision(13) <<second << "       Diverge: " << (second - sqrt(2));
}