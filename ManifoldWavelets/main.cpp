#include "manifoldwavelets.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ManifoldWavelets w;
	w.show();
	return a.exec();
}
