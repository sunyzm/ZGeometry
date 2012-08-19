#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include <QtGui/QMainWindow>
#include "ui_manifoldwavelets.h"
#include "OutputHelper.h"

class QManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:
	friend OutputHelper;
	QManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QManifoldWavelets();

private:
	QString statusBarMsg;
	Ui::ManifoldWaveletsClass ui;
	
	void InfoOutput(const QString& msg, int venue = 1); //1: console, 2: status bar
};

#endif // MANIFOLDWAVELETS_H
