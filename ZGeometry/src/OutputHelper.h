#pragma once
#include <QtGui/QStatusBar>
#include <QtGui/QPlainTextEdit>
#include <string>
#include <vector>

enum OutputVenue {OUT_TERMINAL = 0, OUT_CONSOLE = 1, OUT_STATUS = 2, OUT_MSGBOX = 3};

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
	void output(char c, int repeat, int venue = OUT_CONSOLE);
	void clearOutput(int venue = OUT_CONSOLE);
	void outputDateTime(int venue = OUT_CONSOLE);
private:
	QPlainTextEdit* consoleOutput;
	QStatusBar* statusBar;
};
