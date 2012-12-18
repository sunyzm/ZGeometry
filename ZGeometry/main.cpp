#include "QZGeometry.h"
#include <QtGui/QApplication>
#include <QtGui/QLabel>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	QZGeometry w;
	if (!w.initialize()) exit(-1);
	w.show();
	return a.exec();
}
