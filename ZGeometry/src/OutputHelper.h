#pragma once
#include <string>
#include <vector>
#include <QStatusBar>
#include <QPlainTextEdit>
#include <QLabel>

enum OutputVenue {OUT_TERMINAL = 0, OUT_CONSOLE = 1, OUT_STATUS = 2, OUT_MSGBOX = 3, OUT_STATUS_LABEL = 4};

class OutputHelper
{
public:
	OutputHelper(void);
	~OutputHelper(void);
	void setConsole(QPlainTextEdit* pte) {consoleOutput = pte; }
	void setStatusBar(QStatusBar* sb) { statusBar = sb; }
	void setLabel(QLabel* l) { statusLabel = l; } 

	void output(const QString& msg) { output(msg, mDefaultVenue); }
	void output(const char* msg, int venue, double timeout = 0.0); 	
	void output(const std::string& msg, int venue, double timeout = 0.0); 
    void output(const QString& msg, int venue, double timeout = 0.0);
	void output(char c, int repeat, int venue);
    void outputConsole(const QString& msg) { output(msg, OUT_CONSOLE); }
    void outputStatus(const QString& msg) { output(msg, OUT_STATUS); }
	void clearOutput(int venue);
	void outputDateTime(int venue);
	void setDefaultVenue(int venue) { mDefaultVenue = venue; }

private:
	QPlainTextEdit* consoleOutput;
	QStatusBar* statusBar;
	QLabel* statusLabel;
	int mDefaultVenue;
};
