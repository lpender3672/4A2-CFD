#include "visWidget.h"
#include <iostream>
#include <fstream>
#include <algorithm> // added for min_element
#include <vector>    // added for RAII buffers

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent)
{
    
    mainLayout = new QVBoxLayout(this);
    tabWidget = new QTabWidget(this);

    mainLayout->addWidget(tabWidget);
    connect(tabWidget, &QTabWidget::currentChanged, this, &VisWidget::onTabChanged);

    setupTabs();

    sharedXAxisRange = QCPRange(-1, 1);
    sharedYAxisRange = QCPRange(-1, 1);

}

VisWidget::~VisWidget()
{
}

void VisWidget::setupTabs()
{

    tabNames = {
        "ro",
        "rovx",
        "rovy",
        "roe",
        "T",
        "P",
        "hstag",
        "Mach",
        "dro"
    };
    int tabCount = tabNames.size();

    for (int i = 0; i < tabCount; ++i) {
        QCustomPlot *plot = nullptr;
        QMeshPlot *meshPlot = nullptr;
        QCPColorScale *colorScale = nullptr;

        createMeshGraph(plot, meshPlot, colorScale, tabNames[i], "X", "Y");

        customPlots.append(plot);
        meshPlots.append(meshPlot);
        colorScales.append(colorScale);

        tabWidget->addTab(plot, tabNames[i]);
    }
}

namespace {

void calc_temp(const t_grid *grid, float *temp)
{
    constexpr float cv = 287.05f / (1.4f - 1.0f);
    const int n = grid->ni * grid->nj;
    for (int idx = 0; idx < n; ++idx)
        temp[idx] = (grid->roe[idx] / grid->ro[idx]
                     - 0.5f * (grid->vx[idx]*grid->vx[idx] + grid->vy[idx]*grid->vy[idx])) / cv;
}

void calc_mach(const t_grid *grid, float *mach, const float *temp)
{
    constexpr float gam_rgas = 1.4f * 287.05f;
    const int n = grid->ni * grid->nj;
    for (int idx = 0; idx < n; ++idx)
        mach[idx] = std::sqrt(grid->vx[idx]*grid->vx[idx] + grid->vy[idx]*grid->vy[idx])
                    / std::sqrt(gam_rgas * temp[idx]);
}

} // namespace

void VisWidget::onTabChanged(int index)
{
    if (hasCfState)
        renderCurrentCfTab(index);
    else
        renderCurrentBlockTab(index);
}

void VisWidget::renderCurrentBlockTab(int index)
{
    if (currentGrids.isEmpty()) return;

    sharedXAxisRange = customPlots[lastTabIndex]->xAxis->range();
    sharedYAxisRange = customPlots[lastTabIndex]->yAxis->range();

    QVector<const float *> ptrs;
    std::vector<std::vector<float>> bufs; // scratch for computed fields
    t_data_type dtype = t_data_type::NODE;

    // Collects a direct struct-member pointer from every grid in one line.
    auto addField = [&](float * const t_grid::* field) {
        for (const auto &g : currentGrids) ptrs.push_back(g.*field);
    };

    switch (index) {
    case 0: addField(&t_grid::ro);    break;
    case 1: addField(&t_grid::rovx);  break;
    case 2: addField(&t_grid::rovy);  break;
    case 3: addField(&t_grid::roe);   break;
    case 4:
        for (const auto &g : currentGrids) {
            bufs.emplace_back(g.ni * g.nj);
            calc_temp(&g, bufs.back().data());
            ptrs.push_back(bufs.back().data());
        }
        break;
    case 5: addField(&t_grid::p);     break;
    case 6: addField(&t_grid::hstag); break;
    case 7:
        for (const auto &g : currentGrids) {
            std::vector<float> tmp(g.ni * g.nj), mch(g.ni * g.nj);
            calc_temp(&g, tmp.data());
            calc_mach(&g, mch.data(), tmp.data());
            bufs.push_back(std::move(mch));
            ptrs.push_back(bufs.back().data());
        }
        break;
    case 8:
        addField(&t_grid::dro);
        dtype = t_data_type::CELL;
        break;
    default: return;
    }

    updateBlockMeshGraph(customPlots[index], meshPlots[index], colorScales[index],
                         currentGrids, ptrs, dtype);
    lastTabIndex = index;
}

void VisWidget::createMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, QString title, QString xTitle, QString yTitle)
{
    if (customPlot)
    {
        delete customPlot;
        customPlot = nullptr;
    }

    customPlot = new QCustomPlot(this); // Or in your widget
    meshPlot = new QMeshPlot(customPlot->xAxis, customPlot->yAxis);
    colorScale = new QCPColorScale(customPlot);

    customPlot->plotLayout()->addElement(0, 1, colorScale); // Color scale on the right
    colorScale->setType(QCPAxis::atRight);
    colorScale->setGradient(QCPColorGradient::gpJet);
    customPlot->plotLayout()->setColumnStretchFactor(0, 0.8);
    customPlot->plotLayout()->setColumnStretchFactor(1, 0.2);

    customPlot->xAxis->setLabel(xTitle);
    customPlot->yAxis->setLabel(yTitle);

    customPlot->replot();
}

void VisWidget::finaliseAxes(QCustomPlot *customPlot, double minX, double maxX, double minY, double maxY)
{
    customPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    customPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    if (resetAxis) {
        customPlot->xAxis->setRange(minX, maxX);
        customPlot->yAxis->setRange(minY, maxY);
        customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);
        resetAxis = false;
    } else {
        customPlot->xAxis->setRange(sharedXAxisRange);
        customPlot->yAxis->setRange(sharedYAxisRange);
    }
    customPlot->replot();
}

void VisWidget::drawPolyOverlay(QCustomPlot *customPlot)
{
    if (currentPolyX.isEmpty()) return;

    // Close the polygon by appending the first point
    QVector<double> t(currentPolyX.size() + 1);
    QVector<double> px(currentPolyX.size() + 1);
    QVector<double> py(currentPolyY.size() + 1);
    for (int i = 0; i < currentPolyX.size(); ++i) {
        t[i] = i; px[i] = currentPolyX[i]; py[i] = currentPolyY[i];
    }
    t.back()  = currentPolyX.size();
    px.back() = currentPolyX[0];
    py.back() = currentPolyY[0];

    auto *curve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
    curve->setData(t, px, py);
    curve->setPen(QPen(Qt::white, 1.5));
    curve->setScatterStyle(QCPScatterStyle::ssNone);
}

void VisWidget::updateBlockMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const QVector<t_grid> &grids, const QVector<const float *> &mesh_datas, t_data_type mesh_data_type)
{
    if (grids.isEmpty() || mesh_datas.size() != grids.size()) return;

    customPlot->clearItems();
    meshPlot->clearPolygons();

    // Compute value and spatial ranges over all grids in one pass
    float minVal = std::numeric_limits<float>::max();
    float maxVal = std::numeric_limits<float>::lowest();
    double minX  = std::numeric_limits<double>::max();
    double maxX  = std::numeric_limits<double>::lowest();
    double minY  = std::numeric_limits<double>::max();
    double maxY  = std::numeric_limits<double>::lowest();

    for (int g = 0; g < grids.size(); ++g) {
        const t_grid &grid   = grids[g];
        const float *data    = mesh_datas[g];
        const int nNodes     = grid.ni * grid.nj;
        const int nCells     = (grid.ni - 1) * (grid.nj - 1);
        const int dataCount  = (mesh_data_type == t_data_type::CELL) ? nCells : nNodes;

        minVal = std::min(minVal, *std::min_element(data, data + dataCount));
        maxVal = std::max(maxVal, *std::max_element(data, data + dataCount));
        minX = std::min(minX, (double)*std::min_element(grid.x, grid.x + nNodes));
        maxX = std::max(maxX, (double)*std::max_element(grid.x, grid.x + nNodes));
        minY = std::min(minY, (double)*std::min_element(grid.y, grid.y + nNodes));
        maxY = std::max(maxY, (double)*std::max_element(grid.y, grid.y + nNodes));
    }
    if (minVal >= maxVal) maxVal = minVal + 1.0f;

    QCPRange dataRange(0.0, 1.0);

    for (int g = 0; g < grids.size(); ++g) {
        const t_grid &grid  = grids[g];
        const float *data   = mesh_datas[g];

        for (int i = 0; i < grid.ni - 1; ++i) {
            for (int j = 0; j < grid.nj - 1; ++j) {
                QPolygonF polygon;
                polygon << QPointF(grid.x[ j      * grid.ni + i    ], grid.y[ j      * grid.ni + i    ])
                        << QPointF(grid.x[ j      * grid.ni + i + 1], grid.y[ j      * grid.ni + i + 1])
                        << QPointF(grid.x[(j + 1) * grid.ni + i + 1], grid.y[(j + 1) * grid.ni + i + 1])
                        << QPointF(grid.x[(j + 1) * grid.ni + i    ], grid.y[(j + 1) * grid.ni + i    ]);

                float cellVal;
                if (mesh_data_type == t_data_type::CELL) {
                    cellVal = data[j * (grid.ni - 1) + i];
                } else {
                    cellVal = (data[ j      * grid.ni + i    ] + data[ j      * grid.ni + i + 1] +
                               data[(j + 1) * grid.ni + i    ] + data[(j + 1) * grid.ni + i + 1]) * 0.25f;
                }

                float norm = (cellVal - minVal) / (maxVal - minVal);
                meshPlot->addPolygon(polygon, colorScale->gradient().color(norm, dataRange));
            }
        }
    }

    colorScale->axis()->setRange(minVal, maxVal);
    colorScale->axis()->setLabel(tabNames[tabWidget->currentIndex()]);
    drawPolyOverlay(customPlot);
    finaliseAxes(customPlot, minX, maxX, minY, maxY);
}

void VisWidget::updateLodMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const lod_mesh &mesh)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    // Basic validation
    if (mesh.length <= 0 || mesh.cells == nullptr) {
        // nothing to draw
        return;
    }

    // Compute global bounds and level range
    double globalMinX = std::numeric_limits<double>::max();
    double globalMaxX = std::numeric_limits<double>::lowest();
    double globalMinY = std::numeric_limits<double>::max();
    double globalMaxY = std::numeric_limits<double>::lowest();
    int minLevel = std::numeric_limits<int>::max();
    int maxLevel = std::numeric_limits<int>::lowest();

    for (int i = 0; i < mesh.length; ++i) {
        const auto &c = mesh.cells[i]; // assume mesh.cells is an array of cell2d-like structs
        globalMinX = std::min(globalMinX, (double)c.xmin);
        globalMaxX = std::max(globalMaxX, (double)c.xmax);
        globalMinY = std::min(globalMinY, (double)c.ymin);
        globalMaxY = std::max(globalMaxY, (double)c.ymax);
        minLevel = std::min(minLevel, (int)c.level);
        maxLevel = std::max(maxLevel, (int)c.level);
    }

    QCPRange dataRange(0.0, 1.0);

    // Draw each cell as a rectangle polygon, color by level (if available)
    for (int i = 0; i < mesh.length; ++i) {
        const auto &c = mesh.cells[i];
        QPolygonF polygon;
        polygon << QPointF((double)c.xmin, (double)c.ymin)
                << QPointF((double)c.xmax, (double)c.ymin)
                << QPointF((double)c.xmax, (double)c.ymax)
                << QPointF((double)c.xmin, (double)c.ymax);

        float normalized = 0.5f;
        if (maxLevel > minLevel) {
            normalized = float(c.level - minLevel) / float(maxLevel - minLevel);
        }

        QColor color = Qt::gray;
        if (colorScale) {
            color = colorScale->gradient().color(normalized, dataRange);
        }

        meshPlot->addPolygon(polygon, color);
    }

    // --- New: connect centers in Morton-key order --------------------------
    // Build index array and sort by cell key (Morton)
    
    if (mesh.length > 1) {
    QVector<double> cx(mesh.length), cy(mesh.length), t(mesh.length);
    for (int i = 0; i < mesh.length; ++i) {
        const auto &c = mesh.cells[i];
        cx[i] = 0.5 * (c.xmin + c.xmax);
        cy[i] = 0.5 * (c.ymin + c.ymax);
        t[i]  = i;
    }

    auto *curve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
    curve->setData(t, cx, cy);
    curve->setPen(QPen(Qt::white, 1));
    curve->setName("Hilbert path (sequential)");
    curve->setScatterStyle(QCPScatterStyle::ssNone);
    }
    // -----------------------------------------------------------------------

    if (colorScale) {
        colorScale->axis()->setRange(minLevel, maxLevel);
        colorScale->axis()->setLabel("Level");
    }

    drawPolyOverlay(customPlot);

    customPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    customPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    customPlot->xAxis->setRange(globalMinX, globalMaxX);
    customPlot->yAxis->setRange(globalMinY, globalMaxY);
    customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);

    customPlot->replot();
}

void VisWidget::outputGrid(const t_grid &grid)
{
    currentGrids = { grid };
    onTabChanged(tabWidget->currentIndex());
}

void VisWidget::outputGridVector(const QVector<t_grid> &gridVector)
{
    if (gridVector.isEmpty()) return;
    currentGrids = gridVector;
    onTabChanged(tabWidget->currentIndex());
}

void VisWidget::outputLodMesh(const lod_mesh &mesh)
{
    currentGrids.clear();
    hasCfState = false;

    currentCells.resize(mesh.length);
    for (int i = 0; i < mesh.length; ++i)
        currentCells[i] = mesh.cells[i];

    currentPolyX.resize(mesh.poly_count);
    currentPolyY.resize(mesh.poly_count);
    for (int i = 0; i < mesh.poly_count; ++i) {
        currentPolyX[i] = mesh.poly_x[i];
        currentPolyY[i] = mesh.poly_y[i];
    }

    resetAxis = true;
    updateLodMeshGraph(customPlots[0], meshPlots[0], colorScales[0], mesh);
}

void VisWidget::outputCfState(const cf_state_c &state)
{
    if (state.length <= 0 || state.length != currentCells.size())
        return;

    const int n = state.length;
    cfRo.assign   (state.ro,    state.ro    + n);
    cfRovx.assign (state.rovx,  state.rovx  + n);
    cfRovy.assign (state.rovy,  state.rovy  + n);
    cfRoe.assign  (state.roe,   state.roe   + n);
    cfP.assign    (state.p,     state.p     + n);
    cfVx.assign   (state.vx,    state.vx    + n);
    cfVy.assign   (state.vy,    state.vy    + n);
    cfHstag.assign(state.hstag, state.hstag + n);
    hasCfState = true;

    renderCurrentCfTab(tabWidget->currentIndex());
}

void VisWidget::updateCfMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot,
                                   QCPColorScale *&colorScale,
                                   const QVector<cell2d> &cells,
                                   const std::vector<double> &data,
                                   const QString &label)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    if (cells.isEmpty() || (int)data.size() != cells.size())
        return;

    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();
    for (int i = 0; i < cells.size(); ++i) {
        if (cells[i].iswall) continue;
        minVal = std::min(minVal, data[i]);
        maxVal = std::max(maxVal, data[i]);
    }
    if (minVal >= maxVal) maxVal = minVal + 1.0;

    QCPRange dataRange(0.0, 1.0);

    double globalMinX = std::numeric_limits<double>::max();
    double globalMaxX = std::numeric_limits<double>::lowest();
    double globalMinY = std::numeric_limits<double>::max();
    double globalMaxY = std::numeric_limits<double>::lowest();

    for (int i = 0; i < cells.size(); ++i) {
        const auto &c = cells[i];
        globalMinX = std::min(globalMinX, c.xmin);
        globalMaxX = std::max(globalMaxX, c.xmax);
        globalMinY = std::min(globalMinY, c.ymin);
        globalMaxY = std::max(globalMaxY, c.ymax);

        QPolygonF polygon;
        polygon << QPointF(c.xmin, c.ymin) << QPointF(c.xmax, c.ymin)
                << QPointF(c.xmax, c.ymax) << QPointF(c.xmin, c.ymax);

        QColor color = Qt::gray;
        if (!c.iswall) {
            double norm = (data[i] - minVal) / (maxVal - minVal);
            color = colorScale->gradient().color(norm, dataRange);
        }
        meshPlot->addPolygon(polygon, color);
    }

    colorScale->axis()->setRange(minVal, maxVal);
    colorScale->axis()->setLabel(label);
    drawPolyOverlay(customPlot);
    finaliseAxes(customPlot, globalMinX, globalMaxX, globalMinY, globalMaxY);
}

void VisWidget::renderCurrentCfTab(int index)
{
    if (!hasCfState || currentCells.isEmpty()) return;

    const double gam  = 1.4;
    const double rgas = 287.05;
    const int n = currentCells.size();

    sharedXAxisRange = customPlots[lastTabIndex]->xAxis->range();
    sharedYAxisRange = customPlots[lastTabIndex]->yAxis->range();

    switch (index) {
    case 0: updateCfMeshGraph(customPlots[0], meshPlots[0], colorScales[0], currentCells, cfRo,    "ro");    break;
    case 1: updateCfMeshGraph(customPlots[1], meshPlots[1], colorScales[1], currentCells, cfRovx,  "rovx");  break;
    case 2: updateCfMeshGraph(customPlots[2], meshPlots[2], colorScales[2], currentCells, cfRovy,  "rovy");  break;
    case 3: updateCfMeshGraph(customPlots[3], meshPlots[3], colorScales[3], currentCells, cfRoe,   "roe");   break;
    case 4: {
        std::vector<double> temp(n);
        for (int i = 0; i < n; ++i)
            temp[i] = (cfRo[i] > 0.0) ? cfP[i] / (cfRo[i] * rgas) : 0.0;
        updateCfMeshGraph(customPlots[4], meshPlots[4], colorScales[4], currentCells, temp, "T");
        break;
    }
    case 5: updateCfMeshGraph(customPlots[5], meshPlots[5], colorScales[5], currentCells, cfP,     "p");     break;
    case 6: updateCfMeshGraph(customPlots[6], meshPlots[6], colorScales[6], currentCells, cfHstag, "hstag"); break;
    case 7: {
        std::vector<double> mach(n);
        for (int i = 0; i < n; ++i) {
            double a2 = (cfRo[i] > 0.0) ? gam * cfP[i] / cfRo[i] : 1.0;
            mach[i] = std::sqrt(cfVx[i]*cfVx[i] + cfVy[i]*cfVy[i]) / std::sqrt(std::max(a2, 1e-10));
        }
        updateCfMeshGraph(customPlots[7], meshPlots[7], colorScales[7], currentCells, mach, "Mach");
        break;
    }
    default: break;
    }

    lastTabIndex = index;
}
