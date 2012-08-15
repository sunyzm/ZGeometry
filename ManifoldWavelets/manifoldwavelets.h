#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include <QtGui/QMainWindow>
#include "ui_manifoldwavelets.h"

class ManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:
	ManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~ManifoldWavelets();

private:
	Ui::ManifoldWaveletsClass ui;
};

#endif // MANIFOLDWAVELETS_H
