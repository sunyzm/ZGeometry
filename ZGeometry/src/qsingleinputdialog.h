#ifndef QSINGLEINPUTDIALOG_H
#define QSINGLEINPUTDIALOG_H

#include <QDialog>
#include "ui_qsingleinputdialog.h"

class QSingleInputDialog : public QDialog
{
    Q_OBJECT

public:
    QSingleInputDialog(QWidget *parent = 0);
    ~QSingleInputDialog();

private:
    Ui::QSingleInputDialog ui;
};

#endif // QSINGLEINPUTDIALOG_H
