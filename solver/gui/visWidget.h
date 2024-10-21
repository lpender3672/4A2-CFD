#ifndef VISWIDGET_H
#define VISWIDGET_H

#include <QWidget>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChart>
#include <QVBoxLayout>

#include "../types.h"

class VisWidget : public QWidget
{
    Q_OBJECT

public:

    VisWidget(QWidget *parent = 0);
    ~VisWidget();

    void outputGrid(const t_grid &grid);

signals:
    void newGrid(const t_grid &message);

private:

    QChartView *chartView1;
    QChartView *chartView2;

    void createGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle);
};

#endif // VISWIDGET_H