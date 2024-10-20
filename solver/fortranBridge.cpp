
#include "fortranBridge.h"
#include "gui/consoleWidget.h"

static ConsoleWidget* globalConsoleWidget = nullptr;

void setGlobalConsoleWidget(ConsoleWidget* widget) {
    globalConsoleWidget = widget;

    QObject::connect(globalConsoleWidget, &ConsoleWidget::newMessage,
                     globalConsoleWidget, &ConsoleWidget::outputMessage);
}

extern "C" {

void emit_console_signal(const char* text, int length) {
    if (globalConsoleWidget && length > 0) {
        QString message = QString::fromUtf8(text, length);
        emit globalConsoleWidget->newMessage(message);
    }
}

void emit_grid_signal(t_grid g) {
    if (globalConsoleWidget) {
        emit globalConsoleWidget->newMessage("grid signal received");

        for (int i = 0; i < g.ni; i++) {
            QString row = "";
            for (int j = 0; j < g.nj; j++) {
                row += QString::number(g.x[i * g.ni + j]) + " ";
            }
            emit globalConsoleWidget->newMessage(row);
        }
    }
}
}