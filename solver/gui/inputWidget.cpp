
#include "inputwidget.h"
#include <QFileDialog>
#include <QDebug>

InputWidget::InputWidget(QWidget *parent) : QWidget(parent)
{
    QVBoxLayout *layout = new QVBoxLayout(this);

    choosePathButton = new QPushButton("Choose Path", this);
    runButton = new QPushButton("Run Solver", this);
    pathInput = new QLineEdit(this);
    pathInput->setReadOnly(true);

    layout->addWidget(pathInput);
    layout->addWidget(choosePathButton);
    layout->addWidget(runButton);

    setLayout(layout);

    connect(choosePathButton, &QPushButton::clicked, this, &InputWidget::choosePath);
    connect(runButton, &QPushButton::clicked, this, &InputWidget::runSolverRequested);
}

InputWidget::~InputWidget() {}

QString InputWidget::getPath() const {
    return path;
}

void InputWidget::setPath(const QString &newPath) {
    path = newPath;
    pathInput->setText(newPath);  // Update the input box with the chosen path
}

void InputWidget::choosePath() {
    QString chosenPath = QFileDialog::getOpenFileName(this, "Choose File", "", "All Files (*.*)");
    if (!chosenPath.isEmpty()) {
        setPath(chosenPath);
        emit pathChanged(chosenPath);  // Emit a signal that the path has been chosen
    }
}