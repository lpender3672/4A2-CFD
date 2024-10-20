#include "types.h"
#include <QString>

extern "C" {
    void emit_console_signal(const char* text, int length);
    void emit_grid_signal(t_grid g);
}

void setGlobalConsoleWidget(class ConsoleWidget* widget);
