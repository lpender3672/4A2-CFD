#include "consolewidget.h"

ConsoleWidget::ConsoleWidget(QWidget *parent) : QTextEdit(parent) {
    setReadOnly(true);
}

void ConsoleWidget::outputMessage(const QString &message) {
    this->append(message);
}