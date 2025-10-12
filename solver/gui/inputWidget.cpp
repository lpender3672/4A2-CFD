
#include "inputwidget.h"
#include <QFileDialog>
#include <QDebug>
#include <QGridLayout>
#include <QLabel>
#include "../routines.h"

InputWidget::InputWidget(QWidget *parent) : QWidget(parent)
{
    QVBoxLayout *layout = new QVBoxLayout(this);

    choosePathButton = new QPushButton("Choose Path", this);
    runButton = new QPushButton("Run Solver", this);
    stopButton = new QPushButton("Stop Solver", this);
    pathInput = new QLineEdit(this);
    pathInput->setReadOnly(true);

    stopButton->setEnabled(false);

    layout->addWidget(pathInput);
    layout->addWidget(choosePathButton);
    layout->addWidget(runButton);
    layout->addWidget(stopButton);

    QGridLayout *gridLayout = new QGridLayout();
    QStringList leftLabels = {"CFL", "SFAC", "SFAC_RES", "FACSEC", "FCORR", "NSTEPS"};
    QStringList rightLabels = {"D_MAX", "D_VAR", "ALPHA", "PSTAG", "TSTAG", "POUT"};

    cflInput = new QLineEdit(this);
    sfacInput = new QLineEdit(this);
    sfacResInput = new QLineEdit(this);
    dMaxInput = new QLineEdit(this);
    dVarInput = new QLineEdit(this);
    facSecInput = new QLineEdit(this);
    fCorrInput = new QLineEdit(this);
    nstepsInput = new QLineEdit(this);
    alphaInput = new QLineEdit(this);
    pstagInput = new QLineEdit(this);
    tstagInput = new QLineEdit(this);
    poutInput = new QLineEdit(this);

    QLineEdit *leftInputs[] = {cflInput, sfacInput, sfacResInput, facSecInput, fCorrInput, nstepsInput};
    QLineEdit *rightInputs[] = {dMaxInput, dVarInput, alphaInput, pstagInput, tstagInput, poutInput};

    for (int i = 0; i < leftLabels.size(); ++i) {
        QLabel *label = new QLabel(leftLabels[i], this);
        gridLayout->addWidget(label, i, 0);
        gridLayout->addWidget(leftInputs[i], i, 1);
        // signal
        connect(leftInputs[i], &QLineEdit::editingFinished, this, &InputWidget::saveInputFields);
    }
    for (int i = 0; i < rightLabels.size(); ++i) {
        QLabel *label = new QLabel(rightLabels[i], this);
        gridLayout->addWidget(label, i, 2);
        gridLayout->addWidget(rightInputs[i], i, 3);
        // signal
        connect(rightInputs[i], &QLineEdit::editingFinished, this, &InputWidget::saveInputFields);
    }

    layout->addLayout(gridLayout);

    setLayout(layout);

    connect(choosePathButton, &QPushButton::clicked, this, &InputWidget::choosePath);
    connect(runButton, &QPushButton::clicked, this, &InputWidget::runSolverRequested);
    connect(stopButton, &QPushButton::clicked, this, &InputWidget::stopButtonPressed);
}

InputWidget::~InputWidget() {}

QString InputWidget::getPath() const {
    return path;
}

bool InputWidget::setPath(const QString &newPath) {
    if (newPath.isEmpty()) {
        return false;
    }
    try {
        read_settings(newPath.toStdString(), av, bcs);
    } catch (const std::exception &e) {
        qDebug() << e.what();
        return false;
    }
    updateInputFields();
    path = newPath;
    pathInput->setText(newPath);  // Update the input box with the chosen path
    return true;
}

void InputWidget::choosePath() {
    QString chosenPath = QFileDialog::getOpenFileName(this, "Choose File", "", "All Files (*.*)");
    
    if (setPath(chosenPath)) {
        emit pathChanged(chosenPath);
    }
}

void InputWidget::stopButtonPressed() {
    set_stopit_flag();
}

void InputWidget::updateInputFields() {
    cflInput->setText(QString::number(av.cfl, 'f'));
    sfacInput->setText(QString::number(av.sfac, 'f'));
    sfacResInput->setText(QString::number(av.sfac_res, 'f'));
    dMaxInput->setText(QString::number(av.d_max, 'f'));
    dVarInput->setText(QString::number(av.d_var, 'f'));
    facSecInput->setText(QString::number(av.facsec, 'f'));
    fCorrInput->setText(QString::number(av.fcorr, 'f'));
    nstepsInput->setText(QString::number(av.nsteps));
    alphaInput->setText(QString::number(bcs.alpha, 'f'));
    pstagInput->setText(QString::number(bcs.pstag, 'f'));
    tstagInput->setText(QString::number(bcs.tstag, 'f'));
    poutInput->setText(QString::number(bcs.p_out, 'f'));
}

void InputWidget::blockInputFields() {
    cflInput->setReadOnly(true);
    sfacInput->setReadOnly(true);
    sfacResInput->setReadOnly(true);
    dMaxInput->setReadOnly(true);
    dVarInput->setReadOnly(true);
    facSecInput->setReadOnly(true);
    fCorrInput->setReadOnly(true);
    nstepsInput->setReadOnly(true);
    alphaInput->setReadOnly(true);
    pstagInput->setReadOnly(true);
    tstagInput->setReadOnly(true);
    poutInput->setReadOnly(true);

    choosePathButton->setEnabled(false);
    stopButton->setEnabled(true);
    runButton->setEnabled(false);
}

void InputWidget::unblockInputFields() {
    cflInput->setReadOnly(false);
    sfacInput->setReadOnly(false);
    sfacResInput->setReadOnly(false);
    dMaxInput->setReadOnly(false);
    dVarInput->setReadOnly(false);
    facSecInput->setReadOnly(false);
    fCorrInput->setReadOnly(false);
    nstepsInput->setReadOnly(false);
    alphaInput->setReadOnly(false);
    pstagInput->setReadOnly(false);
    tstagInput->setReadOnly(false);
    poutInput->setReadOnly(false);
    
    choosePathButton->setEnabled(true);
    stopButton->setEnabled(false);
    runButton->setEnabled(true);
}

void InputWidget::saveInputFields() {
    try {
        av.cfl = cflInput->text().toFloat();
        av.sfac = sfacInput->text().toFloat();
        av.sfac_res = sfacResInput->text().toFloat();
        av.d_max = dMaxInput->text().toFloat();
        av.d_var = dVarInput->text().toFloat();
        av.facsec = facSecInput->text().toFloat();
        av.fcorr = fCorrInput->text().toFloat();
        av.nsteps = nstepsInput->text().toInt();
        bcs.alpha = alphaInput->text().toFloat();
        bcs.pstag = pstagInput->text().toFloat();
        bcs.tstag = tstagInput->text().toFloat();
        bcs.p_out = poutInput->text().toFloat();

    } catch (const std::exception &e) {
        qDebug() << e.what();
        return;
    }
    
    write_settings(path.toStdString(), av, bcs);
}
