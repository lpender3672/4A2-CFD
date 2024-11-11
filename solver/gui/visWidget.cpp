#include "visWidget.h"
#include <iostream>

VisWidget::VisWidget(QWidget *parent)
    : QWidget(parent)
{
    
    mainLayout = new QVBoxLayout(this);
    tabWidget = new QTabWidget(this);
    currentGrid = nullptr;

    mainLayout->addWidget(tabWidget);
    connect(tabWidget, &QTabWidget::currentChanged, this, &VisWidget::onTabChanged);

    setupTabs();

}

VisWidget::~VisWidget()
{
    delete currentGrid;
}

void VisWidget::setupTabs()
{
    int tabCount = 8;

    QStringList tabNames = {
        "ro",
        "rovx",
        "rovy",
        "roe",
        "T",
        "P",
        "hstag",
        "mach"
    };

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
    if (currentGrid == nullptr)
    {
        return;
    }

    float *mach = nullptr;
    float *temp = nullptr;

    switch (index)
    {
    case 0:
        updateMeshGraph(customPlots[0], meshPlots[0], colorScales[0], currentGrid, currentGrid->ro, t_data_type::NODE);
        break;
    case 1:
        updateMeshGraph(customPlots[1], meshPlots[1], colorScales[1], currentGrid, currentGrid->rovx, t_data_type::NODE);
        break;
    case 2:
        updateMeshGraph(customPlots[2], meshPlots[2], colorScales[2], currentGrid, currentGrid->rovy, t_data_type::NODE);
        break;
    case 3:
        updateMeshGraph(customPlots[3], meshPlots[3], colorScales[3], currentGrid, currentGrid->roe, t_data_type::NODE);
        break;
    case 4:
        temp = new float[currentGrid->ni * currentGrid->nj];
        calc_temp(currentGrid, temp);
        updateMeshGraph(customPlots[4], meshPlots[4], colorScales[4], currentGrid, temp, t_data_type::NODE);
        break;
    case 5:
        updateMeshGraph(customPlots[5], meshPlots[5], colorScales[5], currentGrid, currentGrid->p, t_data_type::NODE);
        break;
    case 6:
        updateMeshGraph(customPlots[6], meshPlots[6], colorScales[6], currentGrid, currentGrid->hstag, t_data_type::NODE);
        break;
    case 7:
        mach = new float[currentGrid->ni * currentGrid->nj];
        temp = new float[currentGrid->ni * currentGrid->nj];
        calc_temp(currentGrid, temp);
        calc_mach(currentGrid, mach, temp);
        updateMeshGraph(customPlots[7], meshPlots[7], colorScales[7], currentGrid, mach, t_data_type::NODE);
        break;
    default:
        break;
    }
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

void VisWidget::updateMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const t_grid *grid, const float *mesh_data, t_data_type mesh_data_type)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    float minValue, maxValue;

    // the averaged nodes cant be less or greater than the respective min and max node values
    if (mesh_data_type == t_data_type::NODE) {
        minValue = *std::min_element(mesh_data, mesh_data + grid->ni * grid->nj);
        maxValue = *std::max_element(mesh_data, mesh_data + grid->ni * grid->nj);
        // ideally would peform the node averaging here and then pass the averaged data to the meshPlot as cell data
    } else if (mesh_data_type == t_data_type::CELL) {
        minValue = *std::min_element(mesh_data, mesh_data + (grid->ni - 1) * (grid->nj - 1));
        maxValue = *std::max_element(mesh_data, mesh_data + (grid->ni - 1) * (grid->nj - 1));
    } else {
        std::cerr << "Invalid data type" << std::endl;
        return;
    }
    
    float minX = *std::min_element(grid->x, grid->x + grid->ni * grid->nj);
    float maxX = *std::max_element(grid->x, grid->x + grid->ni * grid->nj);
    float minY = *std::min_element(grid->y, grid->y + grid->ni * grid->nj);
    float maxY = *std::max_element(grid->y, grid->y + grid->ni * grid->nj);

    QCPRange dataRange(0, 1);

    for (int i = 0; i < grid->ni - 1; i++)
    {
        for (int j = 0; j < grid->nj - 1; j++)
        {
            QPolygonF polygon;
            float normalizedValue;

            // Define the four corners of the quadrilateral (or triangle if needed)
            polygon << QPointF(grid->x[j * grid->ni + i], grid->y[j * grid->ni + i])        // Bottom-left
                    << QPointF(grid->x[j * grid->ni + i + 1], grid->y[j * grid->ni + i + 1])  // Bottom-right
                    << QPointF(grid->x[(j + 1) * grid->ni + i + 1], grid->y[(j + 1) * grid->ni + i + 1]) // Top-right
                    << QPointF(grid->x[(j + 1) * grid->ni + i], grid->y[(j + 1) * grid->ni + i]);  // Top-left


            if (mesh_data_type == t_data_type::CELL){
                // value is at the center of the cell
                normalizedValue = (mesh_data[j * (grid->ni - 1) + i] - minValue) / (maxValue - minValue);
            } else if (mesh_data_type == t_data_type::NODE){
                // average four corner nodes to get the value at the center of the cell
                normalizedValue = ((mesh_data[j * grid->ni + i] +
                                    mesh_data[j * grid->ni + i + 1] + 
                                    mesh_data[(j + 1) * grid->ni + i] +
                                    mesh_data[(j + 1) * grid->ni + i + 1]) / 4 - minValue) / (maxValue - minValue);
            }

            QColor color = colorScale->gradient().color(normalizedValue, dataRange);
            meshPlot->addPolygon(polygon, color);
        }
    }

    colorScale->axis()->setRange(minValue, maxValue);

    // reset axis
    customPlot->xAxis->setRange(minX, maxX);
    customPlot->yAxis->setRange(minY, maxY);

    // set aspect ratio to 1:1
    customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);

    customPlot->setInteractions(QCP::iRangeZoom | QCP::iRangeDrag);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal | Qt::Vertical);
    customPlot->axisRect()->setRangeDrag(Qt::Horizontal | Qt::Vertical);

    customPlot->replot(); // replot
}

void VisWidget::outputGrid(const t_grid &grid)
{   
    if (currentGrid) {
        *currentGrid = grid; 
    } else {
        currentGrid = new t_grid(grid); 
    }

    onTabChanged(tabWidget->currentIndex());

}
