#ifndef VISWIDGET_H
#define VISWIDGET_H

#include <QWidget>
#include <QVBoxLayout>
#include <QGraphicsLayout>
#include <QTabWidget>
#include <vector>

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
            //painter->setPen(QPen(Qt::black));
            painter->setPen(Qt::NoPen);
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

    explicit VisWidget(QWidget *parent = nullptr);
    ~VisWidget();

    void outputGrid(const t_grid &grid);

    void outputGridVector(const QVector<t_grid> &gridVector);

    void outputLodMesh(const lod_mesh &mesh);

    void outputCfState(const cf_state_c &state);

signals:
    void newGrid(const t_grid &message);

    void newGridVector(const QVector<t_grid> &message);

    void newLodMesh(const lod_mesh &message);

    void newCfState(const cf_state_c &message);

private:
    //StringList tabNames;
    QStringList tabNames;

    QVBoxLayout *mainLayout;
    QTabWidget *tabWidget;

    QVector<QCustomPlot*> customPlots;
    QVector<QCPColorScale*> colorScales;
    QVector<QMeshPlot*> meshPlots;

    QCPRange sharedXAxisRange;
    QCPRange sharedYAxisRange;

    int lastTabIndex = 0;
    bool resetAxis = true;

    void setupTabs();

    QVector<t_grid> currentGrids;

    // CF solver state
    QVector<cell2d> currentCells;
    QVector<double> currentPolyX;
    QVector<double> currentPolyY;
    bool hasCfState = false;
    std::vector<double> cfRo, cfRovx, cfRovy, cfRoe;
    std::vector<double> cfP, cfVx, cfVy, cfHstag;

    void createMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle);
    void finaliseAxes(QCustomPlot *customPlot, double minX, double maxX, double minY, double maxY);
    void drawPolyOverlay(QCustomPlot *customPlot);

    void updateBlockMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const QVector<t_grid> &grids, const QVector<const float *> &mesh_data, t_data_type mesh_data_type);
    void updateLodMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const lod_mesh &mesh);
    void updateCfMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const QVector<cell2d> &cells, const std::vector<double> &data, const QString &label);
    void renderCurrentCfTab(int index);
    void renderCurrentBlockTab(int index);

private slots:
    void onTabChanged(int index);

};

#endif // VISWIDGET_H