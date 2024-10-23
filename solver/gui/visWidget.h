#ifndef VISWIDGET_H
#define VISWIDGET_H

#include <QWidget>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChart>
#include <QVBoxLayout>
#include <QGraphicsLayout>

#include "qcustomplot.h"
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

    QCustomPlot *customPlot1;
    QCustomPlot *customPlot2;

    QCPColorMap *colorMap1;
    QCPColorMap *colorMap2;

    QCPColorScale *colorScale1;
    QCPColorScale *colorScale2;

    t_grid currentGrid;

    void createScatterGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle);
    void createMeshGraph(QCustomPlot *&customPlot, QCPColorMap *&colorMap, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle);

    void updateMeshGraph(QCustomPlot *&customPlot, QCPColorMap *&colorMap, const t_grid &grid, const float *mesh_data);
};

#endif // VISWIDGET_H