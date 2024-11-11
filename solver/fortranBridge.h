#include "types.h"
#include "gui/consoleWidget.h"
#include "gui/visWidget.h"
#include "gui/convWidget.h"
#include <QString>

extern "C" {
    void emit_console_signal(const char* text, int length);
    void emit_grid_signal(t_grid g);
    void emit_conv_point_signal(t_conv_point cp);
}

void setGlobalConsoleWidget(class ConsoleWidget* widget);
void setGlobalVisWidget(class VisWidget* widget);
void setGlobalConvWidget(ConvWidget* widget);
