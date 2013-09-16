#include "zmeshviewer.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ZMeshViewer w;
	w.show();
	return a.exec();
}