
#include "fortranBridge.h"


static ConsoleWidget* globalConsoleWidget = nullptr;
static VisWidget* globalVisWidget = nullptr;

void setGlobalConsoleWidget(ConsoleWidget* widget) {
    globalConsoleWidget = widget;

    QObject::connect(globalConsoleWidget, &ConsoleWidget::newMessage,
                     globalConsoleWidget, &ConsoleWidget::outputMessage);
}

void setGlobalVisWidget(VisWidget* widget) {
    globalVisWidget = widget;

    QObject::connect(globalVisWidget, &VisWidget::newGrid,
                     globalVisWidget, &VisWidget::outputGrid);
}

extern "C" {

// These functions are called by fortran to emit signals to the main GUI thread

void emit_console_signal(const char* text, int length) {
    if (globalConsoleWidget && length > 0) {
        QString message = QString::fromUtf8(text, length);
        emit globalConsoleWidget->newMessage(message);
    }
}

void emit_grid_signal(t_grid g) {

    if (globalVisWidget) {
        emit globalVisWidget->newGrid(g);
    }
    /*
    if (globalConsoleWidget) {
        emit globalConsoleWidget->newMessage("grid signal received");

        for (int i = 0; i < g.ni; i++) {
            QString row = "";
            for (int j = 0; j < g.nj; j++) {
                row += QString::number(g.y[j * g.nj + i]) + " ";
            }
            emit globalConsoleWidget->newMessage(row);
        }

    }
    */
}
}