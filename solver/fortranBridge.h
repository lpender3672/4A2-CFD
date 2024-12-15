
#ifndef FORTRAN_BRIDGE_H
#define FORTRAN_BRIDGE_H

#include "types.h"
#include "gui/consoleWidget.h"
#include "gui/visWidget.h"
#include "gui/convWidget.h"
#include <QString>
#include <atomic>



extern "C" {
    void emit_console_signal(const char* text, int length);
    void emit_grid_signal(t_grid g);
    void emit_conv_point_signal(t_conv_point cp);
    void emit_grid_vector_signal(t_grid *g, int length);

    void set_stopit_flag(bool value);
    extern std::atomic<bool> stopit;
}

void setGlobalConsoleWidget(class ConsoleWidget* widget);
void setGlobalVisWidget(class VisWidget* widget);
void setGlobalConvWidget(class ConvWidget* widget);

#endif // FORTRAN_BRIDGE_H