#pragma once
#include <QtGui/QStatusBar>
#include <QtGui/QPlainTextEdit>
#define OUT2CONSOLE 1
#define OUT2STATUS  2


class OutputHelper
{
public:
	OutputHelper(void);;
	~OutputHelper(void);
	void setConsole(QPlainTextEdit* pte) {consoleOutput = pte; }
	void setStatusBar(QStatusBar* sb) { statusBar = sb; }
	void output(const QString& msg, int venue = 1, double timeout = 0.0); //1: console; 2: status bar
private:
	QPlainTextEdit* consoleOutput;
	QStatusBar* statusBar;
};

