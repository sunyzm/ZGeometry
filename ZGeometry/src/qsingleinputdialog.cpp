#include "qsingleinputdialog.h"

QSingleInputDialog::QSingleInputDialog(QWidget *parent)
    : QDialog(parent)
{
    ui.setupUi(this);
    ui.pushButtonOK->setDefault(true);
}

QSingleInputDialog::~QSingleInputDialog()
{

}
