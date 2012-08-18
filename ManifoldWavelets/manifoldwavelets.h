#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include <QtGui/QMainWindow>
#include "ui_manifoldwavelets.h"

class QManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:
	QManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QManifoldWavelets();

private:
	QString statusBarMsg;
	Ui::ManifoldWaveletsClass ui;
};

#endif // MANIFOLDWAVELETS_H
