#include "manifoldwavelets.h"
#include <QtGui/QApplication>
#include <QtGui/QLabel>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	QManifoldWavelets w;
	w.show();
	return a.exec();
}
