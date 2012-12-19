#pragma once
#include <QtGui/QStatusBar>
#include <QtGui/QPlainTextEdit>
#include <string>
#include <vector>
//#define OUT_CONSOLE 1
//#define OUT_STATUS  2
//#define OUT_MSGBOX  3

enum OutputVenue {OUT_CONSOLE = 1, OUT_STATUS = 2, OUT_MSGBOX = 3};

class OutputHelper
{
public:
	OutputHelper(void);;
	~OutputHelper(void);
	void setConsole(QPlainTextEdit* pte) {consoleOutput = pte; }
	void setStatusBar(QStatusBar* sb) { statusBar = sb; }
	void output(const char* msg, int venue = OUT_CONSOLE, double timeout = 0.0); 
	void output(const QString& msg, int venue = OUT_CONSOLE, double timeout = 0.0); 
	void output(const std::string& msg, int venue = OUT_CONSOLE, double timeout = 0.0); 
	void clearOutput(int venue = OUT_CONSOLE);
	void outputDateTime(int venue = OUT_CONSOLE);
private:
	QPlainTextEdit* consoleOutput;
	QStatusBar* statusBar;
};

