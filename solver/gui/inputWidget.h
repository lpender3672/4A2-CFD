#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include <QWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QVBoxLayout>

class InputWidget : public QWidget
{
    Q_OBJECT

public:
    InputWidget(QWidget *parent = nullptr);
    ~InputWidget();

    QString getPath() const;
    void setPath(const QString &newPath);

signals:
    void pathChanged(const QString &newPath);
    void runSolverRequested();

private slots:
    void choosePath();

private:
    QPushButton *choosePathButton;
    QPushButton *runButton;
    QLineEdit *pathInput;
    QString path;
};

#endif // INPUTWIDGET_H