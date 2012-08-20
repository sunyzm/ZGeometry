#include "OutputHelper.h"
#include <ctime>
#include <QDateTime>
#include <Qtgui/QMessageBox>

OutputHelper::OutputHelper( void ) : consoleOutput(NULL), statusBar(NULL)
{

}


OutputHelper::~OutputHelper(void)
{
}


void OutputHelper::output( const QString& msg, int venue /*= 1*/, double timeout /*= 0.0*/ )
{
	if (venue == OUT_CONSOLE && consoleOutput)
		consoleOutput->insertPlainText(msg + "\n");
	else if (venue == OUT_STATUS && statusBar)
		statusBar->showMessage(msg, timeout);
	else if (venue == OUT_MSGBOX)
	{
		QMessageBox::information(NULL, "Important!", msg, QMessageBox::Ok);
	}
}


void OutputHelper::clearOutput( int venue /*= 1*/ )
{
	if (venue == 1 && consoleOutput)
		consoleOutput->clear();
}

void OutputHelper::outputDateTime( int venue /*= 1*/ )
{
	if (venue == 1 && consoleOutput)
	{
		consoleOutput->insertPlainText("Time stamp: " + QDateTime::currentDateTime().toString("hh:mm:ss  MMM d yyyy\n"));
	}
}