//Pascal Triangle
#include <cmath>
#include <iostream>
unsigned long long fac(unsigned int);
int main() {
	unsigned int N;
	unsigned long long fac1 = 1, fac2 = 1, fac3 = 1;
	unsigned long long coeff = 0;
	std::cout << "Input the number (N) of layers of Pascal Triangle: ";
	std::cin >> N;
	for (unsigned int n = 0; n < N; n++) {
		fac1 *= (n==0 ? 1 : n);
		for (unsigned int i = 0; i <= n; i++) {
			fac2 = fac(i);
			fac3 = fac(n-i);
			coeff = fac1 / fac2 / fac3;
			std::cout << coeff << '\t';
		}
		std::cout << std::endl;
	}
	return 0;
}
unsigned long long fac(unsigned int x){
	unsigned long long f;
	f = ((x==0||x==1) ? 1 : fac(x-1)*x);
	return f;
}