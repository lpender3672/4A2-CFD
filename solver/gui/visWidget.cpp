#include "visWidget.h"
#include <iostream>

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

    float *mach = nullptr;
    float *temp = nullptr;

    QVector<const float *> dataPointers;

    sharedXAxisRange = customPlots[lastTabIndex]->xAxis->range();
    sharedYAxisRange = customPlots[lastTabIndex]->yAxis->range();

    switch (index)
    {
    case 0:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].ro);
        }
        updateMeshGraph(customPlots[0], meshPlots[0], colorScales[0], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 1:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].rovx);
        }
        updateMeshGraph(customPlots[1], meshPlots[1], colorScales[1], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 2:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].rovy);
        }
        updateMeshGraph(customPlots[2], meshPlots[2], colorScales[2], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 3:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].roe);
        }
        updateMeshGraph(customPlots[3], meshPlots[3], colorScales[3], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 4:
        // init size of currentGrids.size null pointers
        for (int i = 0; i < currentGrids.size(); ++i) {
            temp = new float[currentGrids[i].ni * currentGrids[i].nj];
            calc_temp(&currentGrids[i], temp);
            dataPointers.push_back(temp);
        }
        updateMeshGraph(customPlots[4], meshPlots[4], colorScales[4], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 5:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].p);
        }
        updateMeshGraph(customPlots[5], meshPlots[5], colorScales[5], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 6:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].hstag);
        }
        updateMeshGraph(customPlots[6], meshPlots[6], colorScales[6], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 7:
        for (int i = 0; i < currentGrids.size(); ++i) {
            temp = new float[currentGrids[i].ni * currentGrids[i].nj];
            mach = new float[currentGrids[i].ni * currentGrids[i].nj];
            calc_temp(&currentGrids[i], temp);
            calc_mach(&currentGrids[i], mach, temp);
            dataPointers.push_back(mach);
        }
        updateMeshGraph(customPlots[7], meshPlots[7], colorScales[7], currentGrids, dataPointers, t_data_type::NODE);
        break;
    case 8:
        for (int i = 0; i < currentGrids.size(); ++i) {
            dataPointers.push_back(currentGrids[i].dro);
        }
        updateMeshGraph(customPlots[8], meshPlots[8], colorScales[8], currentGrids, dataPointers, t_data_type::CELL);
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

void VisWidget::updateMeshGraph(QCustomPlot *&customPlot, QMeshPlot *&meshPlot, QCPColorScale *&colorScale, const QVector<t_grid> &grids,  const QVector<const float *> &mesh_datas, t_data_type mesh_data_type)
{
    customPlot->clearItems();
    meshPlot->clearPolygons();

    float globalMinValue = std::numeric_limits<float>::max();
    float globalMaxValue = std::numeric_limits<float>::lowest();
    float globalMinX = std::numeric_limits<float>::max();
    float globalMaxX = std::numeric_limits<float>::lowest();
    float globalMinY = std::numeric_limits<float>::max();
    float globalMaxY = std::numeric_limits<float>::lowest();

    for (int g = 0; g < grids.size(); ++g){

        const t_grid &grid = grids[g];
        const float *mesh_data = mesh_datas[g];

        float minValue, maxValue;

        // the averaged nodes cant be less or greater than the respective min and max node values
        if (mesh_data_type == t_data_type::NODE) {
            minValue = *std::min_element(mesh_data, mesh_data + grid.ni * grid.nj);
            maxValue = *std::max_element(mesh_data, mesh_data + grid.ni * grid.nj);
            // ideally would peform the node averaging here and then pass the averaged data to the meshPlot as cell data
        } else if (mesh_data_type == t_data_type::CELL) {
            minValue = *std::min_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
            maxValue = *std::max_element(mesh_data, mesh_data + (grid.ni - 1) * (grid.nj - 1));
        } else {
            std::cerr << "Invalid data type" << std::endl;
            return;
        }
        globalMinValue = std::min(globalMinValue, minValue);
        globalMaxValue = std::max(globalMaxValue, maxValue);
        
        float minX = *std::min_element(grid.x, grid.x + grid.ni * grid.nj);
        float maxX = *std::max_element(grid.x, grid.x + grid.ni * grid.nj);
        float minY = *std::min_element(grid.y, grid.y + grid.ni * grid.nj);
        float maxY = *std::max_element(grid.y, grid.y + grid.ni * grid.nj);

        globalMinX = std::min(globalMinX, minX);
        globalMaxX = std::max(globalMaxX, maxX);
        globalMinY = std::min(globalMinY, minY);
        globalMaxY = std::max(globalMaxY, maxY);
    }

    QCPRange dataRange(0, 1);

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
