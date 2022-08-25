#include "show_progress.hpp"

void show_progress(const long s, const long nsamp, const long n) {
	static const double width = nsamp / (n + 1.0);
	static int next_width = 1;

	if (s == 0) {
		std::cout << "Progress: [";
	}

	if (s > next_width * width) {
		if (next_width % 10 == 0) {
			std::cout << next_width / 10;
		} else {
			std::cout << '.';
		}

		++next_width;
	}

	if (s == nsamp - 1) {
		std::cout << "]\n\n";
	}

	std::cout << std::flush;

	return;
}

/*
    int num_digits(long x) {
	int digits = 0;
	while (x /= 10)
		++digits;
	return digits;
    }

    void show_progress(const long s, const long nsamp) {
	static const int ndigits = num_digits(nsamp);
	const double pc = 100.0 * s / nsamp;

	if (s == nsamp) {
		std::cout << "Progress: 100%\n" << std::endl;
		return;
	}

	std::cout.precision(1);
	std::cout << "Progress: "
		<< std::setw(ndigits)
		<< s + 1 << " / " << nsamp << "   "
		<< std::fixed << std::setw(3)
		<< pc << "\%   "
		<< ((s + 1 < nsamp) ? "\r" : "\n\n")
		<< std::flush;

	// reset flags otherwise final output breaks
	std::cout.precision(6);
	std::cout.unsetf(std::ios_base::floatfield);

	return;
    }
*/


