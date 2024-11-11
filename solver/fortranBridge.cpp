
#include "fortranBridge.h"


static ConsoleWidget* globalConsoleWidget = nullptr;
static VisWidget* globalVisWidget = nullptr;
static ConvWidget* globalConvWidget = nullptr;

void setGlobalConsoleWidget(ConsoleWidget* widget) {
    globalConsoleWidget = widget;

    QObject::connect(globalConsoleWidget, &ConsoleWidget::newMessage,
                     globalConsoleWidget, &ConsoleWidget::outputMessage, Qt::BlockingQueuedConnection);
}

void setGlobalVisWidget(VisWidget* widget) {
    globalVisWidget = widget;

    QObject::connect(globalVisWidget, &VisWidget::newGrid,
                     globalVisWidget, &VisWidget::outputGrid, Qt::BlockingQueuedConnection);
}

void setGlobalConvWidget(ConvWidget* widget) {
    globalConvWidget = widget;

    QObject::connect(globalConvWidget, &ConvWidget::newConvPoint,
                     globalConvWidget, &ConvWidget::outputConvPoint, Qt::BlockingQueuedConnection);
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
}

void emit_conv_point_signal(t_conv_point cp) {
    if (globalConvWidget) {
        emit globalConvWidget->newConvPoint(cp);
    }
}
}