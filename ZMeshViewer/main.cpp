#include "zmeshviewer.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ZMeshViewer w;
	w.show();
	return a.exec();
}
