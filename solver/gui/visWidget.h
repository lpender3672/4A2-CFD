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


// enum for 2d array type cell or node
enum class t_data_type
{
    CELL,
    NODE,
    I_FLUX,
    J_FLUX
};

class QMeshPlot : public QCPAbstractPlottable
{
public:
    QMeshPlot(QCPAxis *keyAxis, QCPAxis *valueAxis)
        : QCPAbstractPlottable(keyAxis, valueAxis)
    {

        int gridSize = 53 * 37;
        polygons.reserve(gridSize);
        colors.reserve(gridSize);
    }

    void addPolygon(const QPolygonF &polygon, const QColor &color)
    {
        polygons.push_back(polygon);
        colors.push_back(color);
    }

    virtual void draw(QCPPainter *painter) override
    {
        // iterate over all polygons and draw them
        for (int i = 0; i < polygons.size(); ++i)
        {
            painter->setPen(QPen(Qt::black));
            painter->setBrush(QBrush(colors[i]));

            // connvert coordinates to pixel coordinates
            QPolygonF pixelPolygon;
            for (const QPointF &point : polygons[i])
            {
                pixelPolygon << coordsToPixels(point.x(), point.y());
            }

            // draw the polygon with pixel coordinates
            painter->drawPolygon(pixelPolygon);
        }
    }

    // abstract methods

    virtual QCPRange getKeyRange(bool &foundRange, QCP::SignDomain inSignDomain) const override
    {
        Q_UNUSED(inSignDomain);
        foundRange = true;
        return QCPRange(0, 1);
    }

    virtual QCPRange getValueRange(bool &foundRange, QCP::SignDomain inSignDomain, const QCPRange &inRange = QCPRange()) const override
    {
        Q_UNUSED(inSignDomain);
        Q_UNUSED(inRange);
        foundRange = true;
        return QCPRange(0, 1);
    }

    virtual void drawLegendIcon(QCPPainter *painter, const QRectF &rect) const override
    {
        Q_UNUSED(painter);
        Q_UNUSED(rect);
    }

    virtual double selectTest(const QPointF &pos, bool onlySelectable, QVariant *details = nullptr) const override
    {
        Q_UNUSED(details);

        if (!onlySelectable || selectable())
        {
            for (const QPolygonF &polygon : polygons)
            {
                if (polygon.containsPoint(pos, Qt::OddEvenFill))
                    return 0.0;  // point is inside a polygon
            }
        }
        return -1.0;  // point is outside all polygons
    }

    void clearPolygons()
    {
        polygons.clear();
        colors.clear();
    }

private:
    std::vector<QPolygonF> polygons;  // Use std::vector for polygons
    std::vector<QColor> colors;       // Use std::vector for colors

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

    QCustomPlot *customPlot1;

    QCPColorScale *colorScale1;

    QMeshPlot *meshPlot1;

    t_grid currentGrid;

    void createScatterGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle);
    void createMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle);

    void updateMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const t_grid &grid, const float *mesh_data, t_data_type type);

};

#endif // VISWIDGET_H