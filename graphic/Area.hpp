#ifndef genomic_Area_h
#define genomic_Area_h

template <typename x_type, typename y_type>
struct Area
{
	x_type xstart, xend;
	y_type ystart, yend;
	
	Area() {}  // no initialization
	
	Area(x_type xb, x_type xe, y_type yb, y_type ye)
	: xstart(xb), xend(xe), ystart(yb), yend(ye) {}
	
	Area& operator+= (const Area& b) {
		xstart += b.xstart;
		xend += b.xend;
		ystart += b.ystart;
		yend += b.yend;
		return *this;
	}
	
	Area& operator-= (const Area& b) {
		xstart -= b.xstart;
		xend -= b.xend;
		ystart -= b.ystart;
		yend -= b.yend;
		return *this;
	}
	
	Area& bound(const Area& b) {
		if (xstart < b.xstart) xstart = b.xstart;
		if (xend > b.xend) xend = b.xend;
		if (ystart < b.ystart) ystart = b.ystart;
		if (yend > b.yend) yend = b.yend;
	}
};

template <typename x_type, typename y_type>
Area<x_type, y_type> operator+(const Area<x_type, y_type>& a, const Area<x_type, y_type>& b) {
	return Area<x_type, y_type>(a) += b;
}

template <typename x_type, typename y_type>
Area<x_type, y_type> operator-(const Area<x_type, y_type>& a, const Area<x_type, y_type>& b) {
	return Area<x_type, y_type>(a) -= b;
}
	
#endif