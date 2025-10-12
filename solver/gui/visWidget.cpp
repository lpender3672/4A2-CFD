#include "visWidget.h"
#include <iostream>
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
    // remove currentGrids
    currentGrids.clear();
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

void calc_temp(const t_grid *grid, float *temp)
{
    float gamma = 1.4;
    float rgas = 287.05;
    float cv = rgas / (gamma - 1);

    // b['t'] = ( b['roe'] / b['ro'] - 0.5 * (b['vx']**2 + b['vy']**2) ) / av['cv']

    for (int i = 0; i < grid->ni; i++) {
        for (int j = 0; j < grid->nj; j++) {
            int idx = j * grid->ni + i;
            temp[idx] = (grid->roe[idx] / grid->ro[idx] - 0.5 * (pow(grid->vx[idx], 2) + pow(grid->vy[idx], 2))) / cv;
        }
    }
}
void calc_mach(const t_grid *grid, float *mach, const float *temp)
{
    float gamma = 1.4;
    float rgas = 287.05;

    // b['mach'] = np.sqrt(b['vx']**2 + b['vy']**2) / np.sqrt( av['gam'] * av['rgas'] * b['t'] )

    for (int i = 0; i < grid->ni; i++) {
        for (int j = 0; j < grid->nj; j++) {
            int idx = j * grid->ni + i;
            mach[idx] = sqrt(pow(grid->vx[idx], 2) + pow(grid->vy[idx], 2)) / sqrt(gamma * rgas * temp[idx]);
        }
    }
}

void VisWidget::onTabChanged(int index)
{

    if (currentGrids.size() == 0)
    {
        return;
    }

    QVector<const float *> dataPointers;

    sharedXAxisRange = customPlots[lastTabIndex]->xAxis->range();
    sharedYAxisRange = customPlots[lastTabIndex]->yAxis->range();

    // RAII buffers to hold computed per-grid arrays while updateBlockMeshGraph runs.
    // Using std::vector ensures memory is freed when this function returns.
    std::vector<std::vector<float>> tempBuffers;
    std::vector<std::vector<float>> machBuffers;

    switch (index)
    {
    case 0:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].ro);
        }
        updateBlockMeshGraph(customPlots[0], meshPlots[0], colorScales[0], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 1:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].rovx);
        }
        updateBlockMeshGraph(customPlots[1], meshPlots[1], colorScales[1], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 2:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].rovy);
        }
        updateBlockMeshGraph(customPlots[2], meshPlots[2], colorScales[2], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 3:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].roe);
        }
        updateBlockMeshGraph(customPlots[3], meshPlots[3], colorScales[3], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 4:
        for (int i = 0; i < currentGrids.size(); ++i) {
            const t_grid &g = currentGrids[i];
            tempBuffers.emplace_back(g.ni * g.nj);
            calc_temp(&g, tempBuffers.back().data());
            dataPointers.push_back(tempBuffers.back().data());
        }
        updateBlockMeshGraph(customPlots[4], meshPlots[4], colorScales[4], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 5:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].p);
        }
        updateBlockMeshGraph(customPlots[5], meshPlots[5], colorScales[5], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 6:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].hstag);
        }
        updateBlockMeshGraph(customPlots[6], meshPlots[6], colorScales[6], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 7:
        for (int i = 0; i < currentGrids.size(); ++i) {
            const t_grid &g = currentGrids[i];
            tempBuffers.emplace_back(g.ni * g.nj);
            machBuffers.emplace_back(g.ni * g.nj);
            calc_temp(&g, tempBuffers.back().data());
            calc_mach(&g, machBuffers.back().data(), tempBuffers.back().data());
            dataPointers.push_back(machBuffers.back().data());
        }
        updateBlockMeshGraph(customPlots[7], meshPlots[7], colorScales[7], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 8:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].dro);
        }
        updateBlockMeshGraph(customPlots[8], meshPlots[8], colorScales[8], currentGrids, dataPointers, t_data_type::CELL);
        break;
    default:
        break;
    }

    lastTabIndex = index;
}

void VisWidget::createScatterGraph(QChartView *chartView, QLineSeries *series, QString title, QString xTitle, QString yTitle)
{
    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle(title);
    chart->createDefaultAxes();
    chart->axisX()->setTitleText(xTitle);
    chart->axisY()->setTitleText(yTitle);

    chart->setMargins(QMargins(0, 0, 0, 0));
    chart->layout()->setContentsMargins(0, 0, 0, 0);
    chart->setBackgroundRoundness(0);

    chartView->setChart(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
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

void VisWidget::updateBlockMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const QVector<t_grid> &grids,  const QVector<const float *> &mesh_datas, t_data_type mesh_data_type)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    // compute global ranges/min/max via helper
    float globalMinValue, globalMaxValue, globalMinX, globalMaxX, globalMinY, globalMaxY;
    if (!computeGlobalRanges(grids, mesh_datas, mesh_data_type,
                             globalMinValue, globalMaxValue,
                             globalMinX, globalMaxX,
                             globalMinY, globalMaxY)) {
        return;
    }

    QCPRange dataRange(0, 1); // normalized data range for color mapping

    for (int g = 0; g < grids.size(); ++g) {
        const t_grid &grid = grids[g];
        const float *mesh_data = mesh_datas[g];

        for (int i = 0; i < grid.ni - 1; i++)
        {
            for (int j = 0; j < grid.nj - 1; j++)
            {
                QPolygonF polygon;
                float normalizedValue;

                // Define the four corners of the quadrilateral (or triangle if needed)
                polygon << QPointF(grid.x[j * grid.ni + i], grid.y[j * grid.ni + i])        // Bottom-left
                        << QPointF(grid.x[j * grid.ni + i + 1], grid.y[j * grid.ni + i + 1])  // Bottom-right
                        << QPointF(grid.x[(j + 1) * grid.ni + i + 1], grid.y[(j + 1) * grid.ni + i + 1]) // Top-right
                        << QPointF(grid.x[(j + 1) * grid.ni + i], grid.y[(j + 1) * grid.ni + i]);  // Top-left


                if (mesh_data_type == t_data_type::CELL){
                    // value is at the center of the cell
                    normalizedValue = (mesh_data[j * (grid.ni - 1) + i] - globalMinValue) / (globalMaxValue - globalMinValue);
                } else if (mesh_data_type == t_data_type::NODE){
                    // average four corner nodes to get the value at the center of the cell
                    normalizedValue = ((mesh_data[j * grid.ni + i] +
                                        mesh_data[j * grid.ni + i + 1] + 
                                        mesh_data[(j + 1) * grid.ni + i] +
                                        mesh_data[(j + 1) * grid.ni + i + 1]) / 4 - globalMinValue) / (globalMaxValue - globalMinValue);
                }

                QColor color = colorScale->gradient().color(normalizedValue, dataRange);
                meshPlot->addPolygon(polygon, color);
            }
        }
    }

    colorScale->axis()->setRange(globalMinValue, globalMaxValue);
    // set colorScale title
    colorScale->axis()->setLabel(
        tabNames[tabWidget->currentIndex()]
    );

    // reset axis
    // set aspect ratio to 1:1

    customPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    customPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    if (resetAxis)
    {
        customPlot->xAxis->setRange(globalMinX, globalMaxX);
        customPlot->yAxis->setRange(globalMinY, globalMaxY);
        customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);
        resetAxis = false;
    }
    else
    {
        customPlot->xAxis->setRange(sharedXAxisRange);
        customPlot->yAxis->setRange(sharedYAxisRange);
    }

    customPlot->replot(); // replot
}

// file-local helper to compute global ranges/min/max for given grids & data
bool VisWidget::computeGlobalRanges(const QVector<t_grid> &grids,
                                    const QVector<const float *> &mesh_datas,
                                    t_data_type mesh_data_type,
                                    float &outGlobalMinValue,
                                    float &outGlobalMaxValue,
                                    float &outGlobalMinX,
                                    float &outGlobalMaxX,
                                    float &outGlobalMinY,
                                    float &outGlobalMaxY)
{
    if (grids.size() == 0 || mesh_datas.size() != grids.size()) {
        std::cerr << "Invalid or empty grids/data" << std::endl;
        return false;
    }

    outGlobalMinValue = std::numeric_limits<float>::max();
    outGlobalMaxValue = std::numeric_limits<float>::lowest();
    outGlobalMinX = std::numeric_limits<float>::max();
    outGlobalMaxX = std::numeric_limits<float>::lowest();
    outGlobalMinY = std::numeric_limits<float>::max();
    outGlobalMaxY = std::numeric_limits<float>::lowest();

    for (int g = 0; g < grids.size(); ++g) {
        const t_grid &grid = grids[g];
        const float *mesh_data = mesh_datas[g];

        float minValue, maxValue;

        if (mesh_data_type == t_data_type::NODE) {
            minValue = *std::min_element(mesh_data, mesh_data + grid.ni * grid.nj);
            maxValue = *std::max_element(mesh_data, mesh_data + grid.ni * grid.nj);
            // ideally would perform node averaging here if needed
        } else if (mesh_data_type == t_data_type::CELL) {
            minValue = *std::min_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
            maxValue = *std::max_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
        } else {
            std::cerr << "Invalid data type" << std::endl;
            return false;
        }

        outGlobalMinValue = std::min(outGlobalMinValue, minValue);
        outGlobalMaxValue = std::max(outGlobalMaxValue, maxValue);

        float minX = *std::min_element(grid.x, grid.x + grid.ni * grid.nj);
        float maxX = *std::max_element(grid.x, grid.x + grid.ni * grid.nj);
        float minY = *std::min_element(grid.y, grid.y + grid.ni * grid.nj);
        float maxY = *std::max_element(grid.y, grid.y + grid.ni * grid.nj);

        outGlobalMinX = std::min(outGlobalMinX, minX);
        outGlobalMaxX = std::max(outGlobalMaxX, maxX);
        outGlobalMinY = std::min(outGlobalMinY, minY);
        outGlobalMaxY = std::max(outGlobalMaxY, maxY);
    }

    return true;
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
    std::vector<int> indices(mesh.length);
    for (int i = 0; i < mesh.length; ++i) indices[i] = i;
    std::sort(indices.begin(), indices.end(),
              [&mesh](int a, int b) {
                  // use unsigned to get sensible ordering for 64-bit keys
                  const unsigned long long ka = static_cast<unsigned long long>(mesh.cells[a].key);
                  const unsigned long long kb = static_cast<unsigned long long>(mesh.cells[b].key);
                  return ka < kb;
              });

    // Compute centers in sorted order
    QVector<double> cx;
    QVector<double> cy;
    cx.reserve(mesh.length);
    cy.reserve(mesh.length);
    for (int idx : indices) {
        const auto &c = mesh.cells[idx];
        cx.append(0.5 * ((double)c.xmin + (double)c.xmax));
        cy.append(0.5 * ((double)c.ymin + (double)c.ymax));
    }

    // Clear previous graphs and add one polyline graph connecting centers
    while (customPlot->graphCount() > 0) {
        customPlot->removeGraph(0);
    }

    if (cx.size() > 1) {
        QCPGraph *g = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);
        g->setData(cx, cy);
        g->setPen(QPen(Qt::black, 1));
        g->setLineStyle(QCPGraph::lsLine);
        g->setScatterStyle(QCPScatterStyle::ssNone);
        g->setName("Morton path");
    }
    // -----------------------------------------------------------------------

    // Configure color scale (use "Level" as the label)
    if (colorScale) {
        colorScale->axis()->setRange(minLevel, maxLevel);
        colorScale->axis()->setLabel("Level");
    }

    // Set axis ranges and aspect ratio
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

    // currentGrids pushed back to be of length 1
    currentGrids.clear();
    currentGrids.push_back(grid);

    onTabChanged(tabWidget->currentIndex());

}

void VisWidget::outputGridVector(const QVector<t_grid> &gridVector)
{
    
    if (gridVector.size() == 0) {
        std::cerr << "Empty grid vector" << std::endl;
        return;
    }

    currentGrids.clear();
    // update currentGrids
    for (int i = 0; i < gridVector.size(); ++i) {
        currentGrids.push_back(gridVector[i]);
    }

    onTabChanged(tabWidget->currentIndex());
}

void VisWidget::outputLodMesh(const lod_mesh &mesh)
{

    currentGrids.clear();

    updateLodMeshGraph(customPlots[0], meshPlots[0], colorScales[0], mesh);

}