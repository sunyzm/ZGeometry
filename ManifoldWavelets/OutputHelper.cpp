#include "OutputHelper.h"


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
