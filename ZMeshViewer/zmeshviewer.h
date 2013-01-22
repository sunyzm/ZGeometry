#ifndef ZMESHVIEWER_H
#define ZMESHVIEWER_H

#include <QtGui/QMainWindow>
#include "ui_zmeshviewer.h"

class ZMeshViewer : public QMainWindow
{
	Q_OBJECT

public:
	ZMeshViewer(QWidget *parent = 0, Qt::WFlags flags = 0);
	~ZMeshViewer();

private:
	Ui::ZMeshViewerClass ui;
};

#endif // ZMESHVIEWER_H
