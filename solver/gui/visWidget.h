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


class MeshPlot : public QCustomPlot
{
public:
    MeshPlot(QWidget *parent = nullptr) : QCustomPlot(parent)
    {
    }

    QVector<QVector<QPointF>> polygons;  // polygons
    QVector<QColor> polygonColors;       // colors for polygons

protected:

    void draw(QCPPainter *painter) override
    {
        QCustomPlot::draw(painter);

        // loop through the polygons and draw black outline and fill with color
        for (int i = 0; i < polygons.size(); ++i)
        {
            QPolygonF qPolygon(polygons[i]); // Convert to QPolygonF for drawing

            // Set the brush to the corresponding color
            painter->setBrush(QBrush(polygonColors[i])); 
            painter->setPen(QPen(Qt::black));  // Polygon outline color

            // Draw the polygon
            painter->drawPolygon(qPolygon);
        }
    }
};



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
    void createMeshGraph(QCustomPlot *&customPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle);

    void updateMeshGraph(QCustomPlot *&customPlot, QCPColorScale *&colorScale, const t_grid &grid, const float *mesh_data);

};

#endif // VISWIDGET_H