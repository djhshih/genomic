#include <cstdio>
#include <cstdlib>

#include "statistics.hpp"

using namespace stats;

int main(int argc, const char* argv[]) {

	if (argc == 1) return 1;

	size_t n = argc-1;
	float* array = new float[n];

	for (size_t i = 0; i < n; ++i) {
		array[i] = std::atof(argv[i+1]);
	}

	std::printf("Array = \n\t");
	for (size_t i = 0; i < n; ++i) {
		std::printf("%f ", array[i]);
	}
	std::printf("\n");

	// descriptive statistics
	std::printf("\nDescriptive statistics\n\n");
	std::printf(
		"Count = %d\n"
		"Mean = %f\n"
		"Variance = \n"
		"\tstandard; two-pass: %f\n"
		"\tone-pass, original: %f\n"
		"\tone-pass, improved: %f\n",
		n,
		mean(array, n),
		var(array, n),
		var2(array, n),
		var3(array, n)
	);

	// datum variance
	std::printf("Datum variance = %f\n", datum_var<float, float*>(array, array+n));

	// running summary
	std::printf("\nRunning statistics\n\n");
	std::printf("\tn\tm\tvar\n");
	summary<float> s;
	for (size_t i = 0; i < n; ++i) {
		s.add(array[i]);
		std::printf("\t%d\t%f\t%f\n",
			s.size(), s.mean(), s.var() );
	}


	// peak detection
	std::printf("\nPeak detection\n\n");
	typedef stats::peak_detector<float, float*> peak_detector;
	peak_detector  detector;
	peak_detector::peaks& p = detector.detect(array, array+n);

	for (size_t i = 0; i < p.size(); i+=2) {
		std::printf("\t%d\t%d\n", p[i], p[i+1]);
	}

	delete array;

	return 0;
}

