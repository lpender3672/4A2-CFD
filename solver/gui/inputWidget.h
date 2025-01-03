#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include <QWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QVBoxLayout>
#include "../types.h"
#include "../fortranBridge.h"

class InputWidget : public QWidget
{
    Q_OBJECT

public:
    InputWidget(QWidget *parent = nullptr);
    ~InputWidget();

    QString getPath() const;
    bool setPath(const QString &newPath);

signals:
    void pathChanged(const QString &newPath);
    void runSolverRequested();

public slots:
    void blockInputFields();
    void unblockInputFields();

private slots:
    void choosePath();
    void stopButtonPressed();

    void updateInputFields();
    void saveInputFields();

private:
    QPushButton *choosePathButton;
    QPushButton *runButton;
    QPushButton *stopButton;
    QLineEdit *pathInput;
    QString path;

    t_appvars av;
    t_bconds bcs;

    QLineEdit *cflInput;
    QLineEdit *sfacInput;
    QLineEdit *sfacResInput;
    QLineEdit *dMaxInput;
    QLineEdit *dVarInput;
    QLineEdit *facSecInput;
    QLineEdit *fCorrInput;
    QLineEdit *nstepsInput;

    QLineEdit *alphaInput;
    QLineEdit *pstagInput;
    QLineEdit *tstagInput;
    QLineEdit *poutInput;
};

#endif // INPUTWIDGET_H