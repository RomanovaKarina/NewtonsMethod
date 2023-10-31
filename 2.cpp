#include <iostream>;
#include <iomanip>;
#include<algorithm>
#include <math.h>
using namespace std;

const double e1 = pow(10, -9);
const double e2 = pow(10, -9);
const int n = 2;


double func1(double x1, double x2) {
	return cos(0.4 * x2 + x1 * x1) + x2 * x2 + x1 * x1 - 1.6;
}
double func2(double x1, double x2) {
	return 1.5 * x1 * x1 - ((x2 * x2) / 0.36) - 1;
}
double dFunc1x1(double x1, double x2) {
	return -2 * x1 * sin(0.4 * x2 + x1 * x1) + 2 * x1;
}
double dFunc1x2(double x1, double x2) {
	return -0.4 * sin(0.4 * x2 + x1 * x1) + 2 * x2;
}
double dFunc2x1(double x1, double x2) {
	return 3 * x1;
}
double dFunc2x2(double x1, double x2) {
	return -5.555 * x2;
}

void Iteration_Newton(double x1, double x2);
void solutionM(double x1, double x2, double M);
double* methodGauss(double mA[n][n], double*);
double norma(double vF);


int main() {
	double x1_1 = 1, x2_1 = -1, x1_2 = -1, x2_2 = 1;

	Iteration_Newton(x1_1, x2_1);
	Iteration_Newton(x1_2, x2_2);
	solutionM(x1_1, x2_1, 0.01);
	solutionM(x1_2, x2_2, 0.01);
	solutionM(x1_1, x2_1, 0.05);
	solutionM(x1_1, x2_1, 0.05);
	solutionM(x1_1, x2_1, 0.1);
	solutionM(x1_1, x2_1, 0.1);
}

void Iteration_Newton(double x1, double x2) {

	int k = 1;
	cout << "k" << setw(12) << "l1" << setw(14) << "l2" << setw(16) << "x1" << setw(14) << "x2" << endl;
	cout << "------------------------------------------------------------------" << endl;

	const int NIT = 5;
	double* dXk = new double[n];
	double dXk1 = 0, dXk2 = 0;
	double l1 = 1, l2 = 1, l1_2 = 0, l2_2 = 0;
	cout << k << "\t " << l1 << " \t" << l2 << " \t" << x1 << " \t" << x2 << endl;
	try{
		while (l1 > e1 && l2 > e2) {
			double Fx[] = { -func1(1,-1), -func2(-1,1) };
			double Jacobian[n][n] = {
				{dFunc1x1(1,-1), dFunc1x2(1,-1)},
				{dFunc2x1(1,-1), dFunc2x2(1,-1)}
			};
			dXk = methodGauss(Jacobian, Fx);
			dXk1 = dXk[0] + x1;
			dXk2 = dXk[1] + x2;
			l1 = norma(func1(x1, x2));
			l1_2 = norma(func2(x1, x2));
			if (l1_2 > l1) {
				l1 = l1_2;
			}
			dXk1 >= 1 ?
				l2 = norma(dXk1 - x1) / dXk1 :
				l2 = norma(dXk1 - x1);
			dXk2 >= 1 ?
				l2_2 = norma(dXk1 - x2) / dXk1 :
				l2_2 = norma(dXk1 - x2);
			if (l2_2 > l2)
				l2 = l2_2;
			x1 = dXk1; x2 = dXk2;
			dXk[0] = x1; dXk[1] = x2;
			cout << k << "\t " << l1 << " \t" << l2 << " \t" << x1 << " \t" << x2 << endl;
			k++;
			if (k > NIT)
			{

				throw 1;
				
			}
		}
	}
	catch (int err) {
		cout << endl;
		cout << "IER=2" << endl;

	}
	k = 0;

};
double* methodGauss(double mA[n][n], double cB[n]) {
	double X[n] = { 0 };
	for (int i = 0; i < n; i++) {
		int maxIndex = i;
		double max = mA[i][i];
		for (int j = i + 1; j < n; j++) {
			if (abs(max) < abs(mA[j][i])) {
				maxIndex = j;
				max = mA[j][i];
			}
		}
		if (i != maxIndex) {
			double root = cB[i];
			cB[i] = cB[maxIndex];
			cB[maxIndex] = root;
			for (int j = 0; j < n; j++) {
				double x = mA[i][j];
				mA[i][j] = mA[maxIndex][j];
				mA[maxIndex][j] = x;
			}
		}
		double a = mA[i][i];
		for (int j = i; j < n; j++)
		{
			mA[i][j] /= a;
		}
		cB[i] /= a;
		for (int j = i + 1; j < n; j++)
		{
			double s = mA[j][i];
			for (int k = i; k < n; k++)
			{
				mA[j][k] -= s * mA[i][k];
			}
			cB[j] -= s * cB[i];
		}
	}

	for (int k = n - 1; k >= 0; k--)
	{
		X[k] = cB[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= mA[k][i] * X[i];
		}
	}
	cout << endl;
	return X;
}

double norma(double vF) {
	double normF = abs(vF);
	normF = max(vF, normF);

	return normF;
}

void solutionM(double x1, double x2, double M) {
	int k = 1;
	cout << "k" << setw(12) << "l1" << setw(14) << "l2" << setw(16) << "x1" << setw(14) << "x2" << setw(14) << "M" << endl;
	cout << "------------------------------------------------------------------" << endl;

	const int NIT = 10;
	double* dXk = new double[n];
	double dXk1 = 0, dXk2 = 0;
	double l1 = 1, l2 = 1, l1_2 = 0, l2_2 = 0;
	try {
		while (l1 > e1 && l2 > e2) {
			double Fx[] = { -func1(1,-1), -func2(-1,1) };
			double Jacobian[n][n] = {
				{(func1(x1 + M * x1, x2) - func1(x1, x2)) / (M * x1), (func1(x1, x2 + M * x2) - func1(x1, x2)) / (M * x2)},
				{(func2(x1 + M * x1, x2) - func2(x1, x2)) / (M * x1), (func2(x1, x2 + M * x2) - func2(x1, x2)) / (M * x2)}
			};
			dXk = methodGauss(Jacobian, Fx);
			dXk1 = dXk[0] + x1;
			dXk2 = dXk[1] + x2;
			l1 = norma(func1(x1, x2));
			l1_2 = norma(func2(x1, x2));
			if (l1_2 > l1) {
				l1 = l1_2;
			}
			dXk1 >= 1 ?
				l2 = norma(dXk1 - x1) / dXk1 :
				l2 = norma(dXk1 - x1);
			dXk2 >= 1 ?
				l2_2 = norma(dXk1 - x2) / dXk1 :
				l2_2 = norma(dXk1 - x2);
			if (l2_2 > l2)
				l2 = l2_2;
			x1 = dXk1; x2 = dXk2;
			dXk[0] = x1; dXk[1] = x2;
			cout << k << "\t " << l1 << " \t" << l2 << " \t" << x1 << " \t" << x2 << " \t" << M << endl;
			k++;
			if (k >= NIT)
			{
				throw 1;
			}
		}
	}
	catch(int e) {
		cout << endl;
		cout << "IER=2" << endl;
	}
	k = 0;

}
//1 -1 -1 1