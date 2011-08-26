#ifndef genomic_Tween_h
#define genomic_Tween_h

namespace tween
{
	// t = current time
	// b = start value
	// c = change in value
	// d = duration
	
	
	template<typename V, typename T> inline
	V linear(T t, V b, V c, T d) {
		return c*t/d + b;
	}
	
	template<typename V, typename T> inline
	V easeQuad(T t, V b, V c, T d) {
		t /= d/2;
		if (t < 1) return c/2 * t*t + b;
		t--;
		return -c/2 * (t*(t-2) - 1) + b;
	}
	
	template<typename V, typename T> inline
	V easeInQuad(T t, V b, V c, T d) {
		t /= d;
		return c*t*t + b;
	};
	
	template<typename V, typename T> inline
	V easeOutQuad(T t, V b, V c, T d) {
		t /= d;
		return -c * t*t(1-2) + b;
	}
	
	
	template<typename V, typename T> inline
	V easeCubic(T t, V b, V c, T d) {
		t /= d/2;
		if (t < 1) return c/2 * t*t*t + b;
		t -= 2;
		return c/2*(t*t*t + 2) + b;
	}
	
	template<typename V, typename T> inline
	V easeInCubic(T t, V b, V c, T d) {
		t /= d;
		return c*t*t*t + b;
	}
	
	template<typename V, typename T> inline
	V easeOutCubic(T t, V b, V c, T d) {
		t /= d;
		t--;
		return c*(t*t*t + 1) + b;
	}
	
	
	/*
	template<typename V, typename T> inline
	V ease(T t, V b, V c, T d) {
	}
	template<typename V, typename T> inline
	V ease(T t, V b, V c, T d) {
	}
	template<typename V, typename T> inline
	V ease(T t, V b, V c, T d) {
	}
	template<typename V, typename T> inline
	V ease(T t, V b, V c, T d) {
	}
	*/
	
}

#endif
