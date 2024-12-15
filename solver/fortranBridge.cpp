
#include "fortranBridge.h"
#include <iostream>
#include <atomic>

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

    QObject::connect(globalVisWidget, &VisWidget::newGridVector,
                        globalVisWidget, &VisWidget::outputGridVector, Qt::BlockingQueuedConnection);
}

void setGlobalConvWidget(ConvWidget* widget) {
    globalConvWidget = widget;

    QObject::connect(globalConvWidget, &ConvWidget::newConvPoint,
                     globalConvWidget, &ConvWidget::outputConvPoint, Qt::BlockingQueuedConnection);
}

extern "C" {

static std::atomic<bool> stopit(false);

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

void emit_grid_vector_signal(t_grid *g, int length) {
    if (globalVisWidget) {
        QVector <t_grid> gridVector;
        for (int n = 0; n < length; n++) {
            gridVector.push_back(g[n]);

            /*
            int i = g[n].ni - 1;
            for (int j = 0; j < g[n].nj; j++) {
                int idx = j * g[n].ni + i;
                std::cout << "x[" << idx << "] = " << g[n].x[idx] << std::endl;
            }
            */
        }
        emit globalVisWidget->newGridVector(gridVector);
    }
}

void emit_conv_point_signal(t_conv_point cp) {
    if (globalConvWidget) {
        emit globalConvWidget->newConvPoint(cp);
    }
}

void set_stopit_flag(bool value) {
    stopit.store(value, std::memory_order_relaxed); // Update the flag
    if (value) {
        qDebug() << "Stopit flag set to true.";
    }
}

}