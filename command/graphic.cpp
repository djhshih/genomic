#include "graphic.hpp"

Window window;

void jumpTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
	window.graph->jumpTo(xstartDest, xendDest, ystartDest, yendDest);
}

void moveTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
	window.graph->moveTo(xstartDest, xendDest, ystartDest, yendDest);
}

void easeTo(x_type xstartDest, x_type xendDest, y_type ystartDest, y_type yendDest) {
	window.graph->easeTo(xstartDest, xendDest, ystartDest, yendDest);
}
