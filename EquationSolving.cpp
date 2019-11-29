#include <iostream>
#include <iomanip>
#include <discpp.h>
#include <stdio.h>

namespace EqSol {
	Dislin g;
	const double eps = 1e-14;
	const double dx = 1e-7;
	double derivative(double (*func)(double), double x)
	{

		return (func(x + dx) - func(x - dx)) / (2*dx);
	}
	double* derivative(const double *y, const double *x, const int n)
	{
		double t=0;
		double *res = new double[n-1];
		for (int i = 0; i < n-1; i++)
			res[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		return res;
	}
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
	double Newton(const double x)
	{
		if (f(x) < eps)
			return x;
		double res = x;
		while (f(res) > eps)
			res = res - f(res) / derivative(f, res);
		return res;
	}
	void secondderiv()
	{	
		g.metafl("cons");
		g.scrmod("revers");
		g.disini();
		g.titlin("Second derivative", 2);
		g.title();
		const int n = 200;
		double *x = new double[n];
		double *t = new double[n];
		double *y = new double[n - 1];
		const double dx = 4.0 / n;
		for (int i = 0; i < n; i++)
		{
			x[i] = -2 + dx * i;
			t[i] = derivative(f, x[i]);
		}
		y = derivative(t, x, n);
		g.qplot(x, y, n - 1);
		remove("dislin.met");
	}
}
int main()
{
	//std::cout << EqSol::derivative(EqSol::f, 10+1E-14);
	double firstb = EqSol::binsolution(-2, 0);
	double secondb = EqSol::binsolution(0, 2);
	double firstn = EqSol::Newton(-2);
	double secondn = EqSol::Newton(2);
	std::cout << "-------------------------------------------------------------------" << '\n';
	std::cout << "Binary solution:" << '\n';
	std::cout << std::setprecision(13) <<firstb <<"      Diverge: "<<abs(firstb+sqrt(2))<<'\n';
	std::cout << std::setprecision(13) <<secondb << "       Diverge: " << (secondb - sqrt(2))<<'\n';
	std::cout << "Newton solution:" << '\n';
	std::cout << std::setprecision(13) <<firstn << "      Diverge: " << abs(firstn + sqrt(2)) << '\n';
	std::cout << std::setprecision(13) << secondn << "       Diverge: " << abs(secondn - sqrt(2)) << '\n';
	std::cout << "-------------------------------------------------------------------" << '\n';
	std::cout << "Step in calculating of first derivation is set at " << EqSol::dx << '\n';
	std::cout << "-------------------------------------------------------------------" << '\n';
	EqSol::secondderiv();
	return 0;
}