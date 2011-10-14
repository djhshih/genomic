#ifndef _statistics_h_
#define _statistics_h_

#include <cstddef>
#include <cmath>
#include <vector>

#include <cstdio>

#include <boost/math/distributions/normal.hpp>

namespace stats {

// arithmetic mean
template <typename T>
T mean(const T* array, size_t n) {
	T r = 0;
	for (size_t i = 0; i < n; ++i) {
		r += array[i];
	}
	return r/n;
}

// standard deviation
// standard (two-pass) method
// one-pass to calculate mean, and another to calcuate var
template <typename T>
T var(const T* array, size_t n, bool sample = true) {
	if (n <= 1) return 0;

	T m = mean(array, n);
	T v = 0;
	for (size_t i = 0; i < n; ++i) {
		T d = array[i] - m;
		v += d*d;
	}
	return sample ? v/(n-1) : v/n;
}

// standard deviation
// one-pass method, original (sum of squares)
template <typename T>
T var2(const T* array, size_t n, bool sample = true) {
	if (n <= 1) return 0;

	T sum = 0;
	T sqsum = 0;
	for (size_t i = 0; i < n; ++i) {
		sum += array[i];
		sqsum += array[i] * array[i];
	}
	return sample
		? sqsum/(n-1) - (sum/(n-1) * sum/n)
		: sqsum/n - (sum/n * sum/n);
}

// standard deviation
// one-pass method, improved (Welford, 1962)
template <typename T>
T var3(const T* array, size_t n, bool sample = true) {
	if (n <= 1) return 0;

	T m = array[0];
	T q = 0;
	for (size_t i = 1; i < n; ++i ) {
		//T t = (array[i] - m);
		//q += i * t*t / (i+1);
		//m += (array[i] - m) / (i+1);

		T t = (array[i] - m);
		m += (array[i] - m) / (i+1);
		q += t * (array[i] - m);
	}

	return sample ? q/(n-1) : q/n;
}

/*

	 For 1 < x < 10, one-pass original gives accurate numbers
	 Beyond this range, if x is float, one-pass original can become 
	   very inaccurate. If x is double, one-pass original is still accurate.
	 
	 Example: 10001, 10002, 10003

*/


// Running statistics
// calculates SD based on Welford method
template <typename T>
class summary
{
public:

	summary() : m(0), devsq(0), n(0) {}

	void add(T x) {
		T t = (x - m);
		m += (x - m) / (++n);
		devsq += t * (x - m);
	}

	size_t size() {
		return n;
	}

	// mean
	T mean() {
		return m;
	}

	// sample standard deviation
	T sd() {
		return (n > 1) ? std::sqrt(devsq/(n-1)) : 0;
	}

	// sample variance
	T var() {
		return (n > 1) ? devsq/(n-1) : 0;
	}

	// population variance
	T pvar() {
		return (n > 1) ? devsq/(n-1) : 0;
	}

	// population standard deviation
	T psd() {
		return (n > 1) ? std::sqrt(devsq/n) : 0;
	}

private:
	T m, devsq;
	size_t n;
};

// datum variance estimation
// assume all data points have the same true variance
// estimate this variance from differences between adjacent data points
template <typename T, typename iterator>
T datum_var(iterator begin, iterator end) {
	summary<T> s;
	for (iterator curr = begin+1; curr != end; ++curr) {
		// add data point difference to running summary statistics
		s.add( (*curr) - (*(curr-1)) );
	}
	// it can be shown that:
	// datum variance to simply the variance of adjacent differences divided by 2
	return s.var() / 2;
}

// peak detector

template < typename T = float,
	         typename iterator = typename std::vector<T>::const_iterator >
class peak_dector;

template <typename T, typename iterator>
class peak_detector
{
public:

	typedef std::vector<size_t> peaks;

	peaks& detect(iterator begin, iterator end, T alpha = 0.05, T region_conf = 0.95) {

		_peaks.clear();

		using namespace boost::math;

		//// determine datum variance, sample variance, and sample mean
		////   in the first pass through the data

		summary<T> diff_summary, data_summary;
		// add the first data point
		data_summary.add(*begin);
		for (iterator curr = begin+1; curr != end; ++curr) {
			// add data point difference to running summary statistics
			diff_summary.add( (*curr) - (*(curr-1)) );
			// add data point to running summary statistics
			data_summary.add(*curr);
		}

		T mean = data_summary.mean();
		T variance = data_summary.var();

		// it can be shown that:
		// datum variance to simply the variance of adjacent differences divided by 2
		T dvariance =  diff_summary.var() / 2;

		//// use sample variance and mean to determine the threshold
		
		// z distribution
		normal_distribution<float> s;

		// critical z-value
		T z_critical = quantile(s, 1 - alpha);

		// threshold
		T threshold = mean + ( z_critical * std::sqrt(variance / data_summary.size()) );

		// reduce threshold to account for datum error
		// (reduce threshold by z * (datum variance)
		// use two-tailed confidence interval
		// for instance, if region_conf == 0.95, convert region_conf to 0.975
		//threshold -= quantile(s, 0.5 + region_conf/2 ) * dvariance;

		
		std::printf("mean = %f\n"
				        "variance = %f\n"
			          "datum error = %f\n"
			          "threshold = %f\n",
			          mean, variance, dvariance, threshold);

		//// report significant regions, in a second pass
		
		size_t i = 0;
		bool in_peak = false;
		for (iterator curr = begin; curr != end; ++curr) {
			if (in_peak) {
				// look for end of peak 	
				if (*curr <= threshold) {
					// current position not in peak any more
					// the previous position is the last position in peak
					_peaks.push_back(i-1);
					in_peak = false;
				}
			} else {
				// look for beginning of peak
				if (*curr > threshold) {
					_peaks.push_back(i);
					in_peak = true;
				}
			}
			++i;
		}

		if (in_peak) {
			// final peak has not been resolved
			// treat last position as end of peak
			_peaks.push_back(i-1);
		}

		return _peaks;
		
	}

private:

	// start and end indices for peaks (all together)
	peaks _peaks;

};


}

#endif
