#include <QTextEdit>

#ifndef CONSOLEWIDGET_H
#define CONSOLEWIDGET_H

class ConsoleWidget : public QTextEdit
{
    Q_OBJECT
public:
    ConsoleWidget(QWidget *parent = nullptr);
    void outputMessage(const QString &message);

signals:
    void newMessage(const QString &message);
};

#endif