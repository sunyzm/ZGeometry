#include "OutputHelper.h"
#include <ctime>
#include <QDateTime>

OutputHelper::OutputHelper( void ) : consoleOutput(NULL), statusBar(NULL)
{

}


OutputHelper::~OutputHelper(void)
{
}


void OutputHelper::output( const QString& msg, int venue /*= 1*/, double timeout /*= 0.0*/ )
{
	if (venue == 1 && consoleOutput)
		consoleOutput->insertPlainText(msg + "\n");
	else if (venue == 2 && statusBar)
		statusBar->showMessage(msg, timeout);
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
// 		time_t curTime;
// 		time(&curTime);
// 		QDateTime qdt;
// 		qdt.setTime_t((uint)curTime);
		consoleOutput->insertPlainText("Time stamp: " + QDateTime::currentDateTime().toString("hh:mm:ss  MMM d yyyy"));
	}
}