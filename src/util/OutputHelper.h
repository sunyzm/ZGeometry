#pragma once
#include <QtGui/QStatusBar>
#include <QtGui/QPlainTextEdit>
#define OUT_CONSOLE 1
#define OUT_STATUS  2


class OutputHelper
{
public:
	OutputHelper(void);;
	~OutputHelper(void);
	void setConsole(QPlainTextEdit* pte) {consoleOutput = pte; }
	void setStatusBar(QStatusBar* sb) { statusBar = sb; }
	void output(const QString& msg, int venue = 1, double timeout = 0.0); //1: console; 2: status bar
	void clearOutput(int venue = 1);
	void outputDateTime(int venue = 1);
private:
	QPlainTextEdit* consoleOutput;
	QStatusBar* statusBar;
};

