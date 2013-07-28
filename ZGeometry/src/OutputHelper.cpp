#include "OutputHelper.h"
#include <ctime>
#include <QDateTime>
#include <Qtgui/QMessageBox>
#include <QScrollBar>
#include <QDebug>
#include <fstream>
#include <iostream>

OutputHelper::OutputHelper( void ) : consoleOutput(NULL), statusBar(NULL) {}

OutputHelper::~OutputHelper(void) {}

void OutputHelper::output( const QString& msg, int venue /*= 1*/, double timeout /*= 0.0*/ )
{
	switch (venue)
	{
		case OUT_TERMINAL: 
			std::cout << msg.toStdString() << std::endl;
//			qDebug() << msg;
			break;
		case OUT_CONSOLE:
			{
				consoleOutput->insertPlainText(msg + "\n");
				QScrollBar *vScrollBar = consoleOutput->verticalScrollBar();
				vScrollBar->triggerAction(QScrollBar::SliderToMaximum);
				break;
			}
		case OUT_STATUS:
			if (statusBar) 	statusBar->showMessage(msg, timeout);
			break;
		case OUT_MSGBOX:
			QMessageBox::information(NULL, "Important!", msg, QMessageBox::Ok);
			break;
		default:
			std::cout << msg.toStdString() << std::endl;
	}	
}

void OutputHelper::output( const std::string& msg, int venue /*= OUT_CONSOLE*/, double timeout /*= 0.0*/ )
{
	output(QString(msg.c_str()), venue, timeout);
}

void OutputHelper::output( const char* msg, int venue /*= OUT_CONSOLE*/, double timeout /*= 0.0*/ )
{
	output(QString(msg), venue, timeout);
}

void OutputHelper::output( char c, int repeat, int venue )
{
	char *str = new char[repeat+1];
	for (int i = 0; i < repeat; ++i) str[i] = c;
	str[repeat] = '\0';
	output(QString(str), venue);
	delete []str;
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

