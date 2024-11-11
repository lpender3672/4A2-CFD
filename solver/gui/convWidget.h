#ifndef CONVWIDGET_H
#define CONVWIDGET_H

#include <QWidget>
#include "../types.h"
#include "qcustomplot.h"



class ConvWidget : public QWidget
{
    Q_OBJECT

public:

    explicit ConvWidget(QWidget *parent = nullptr);
    ~ConvWidget();

    void outputConvPoint(const t_conv_point &cp);

signals:
    void newConvPoint(const t_conv_point &message);

private:
    QVBoxLayout *mainLayout; // No error should occur here now

    QCustomPlot *convPlot;

    QVector<int> iterations;
    QVector<float> d_max;
    QVector<float> d_avg;

    void updateConvGraph();

};

#endif // VISWIDGET_H