from bokeh.io import show, output_file
from bokeh.models import (Circle, MultiLine,
                          GraphRenderer, StaticLayoutProvider, NodesAndLinkedEdges,
                          Div, Column, Row,
                          HoverTool, TapTool, ZoomInTool, ZoomOutTool, PanTool,
                          CheckboxButtonGroup, CustomJS)
from bokeh.layouts import column
from bokeh.plotting import figure, curdoc
from bokeh.colors import RGB
from bokeh.palettes import Spectral4
from bokeh.events import Tap
import parseTests
import testRunner

def makeNetwork():
    nExpRef      = data.nSources
    nTestSeries  = data.nSeries
    nSamples     = data.nSamples
    totalNodes   = nExpRef + nTestSeries + nSamples
    node_indices = list(range(totalNodes))

    x = [0 for i in node_indices]
    y = [0 for i in node_indices]
    names = ['' for i in node_indices]
    types = ['' for i in node_indices]
    details = ['' for i in node_indices]
    fullDetails = ['' for i in node_indices]
    size = [10 for i in node_indices]
    colorCode = [RGB(100,100,100) for i in node_indices]

    level1   = 0.75
    level2   = 0
    level3   = -0.75
    leftLim  = -0.95
    rightLim = 0.95

    compString = 'composition'
    ppString   = 'partial pressures'
    frString   = 'fractions'
    nlString   = '\n'

    edgeStarts = []
    edgeEnds   = []
    sourceIndex = -1
    seriesIndex = nExpRef - 1
    sampleIndex = nExpRef + nTestSeries - 1
    samplePos   = leftLim
    sampleSpace = 0.01
    for sourceName in data.data['sources']:
        source = data.data['sources'][sourceName]
        sourceIndex += 1
        # Calculate node reference positions
        y[sourceIndex] = level1
        x[sourceIndex] = 0
        nExpInSource = 0
        colorCode[sourceIndex] = RGB(0,0,255)
        size[sourceIndex] = 20
        names[sourceIndex] = source['title']
        types[sourceIndex] = source['type']
        details[sourceIndex] = f"Authors: {source['authors']}"
        fullDetails[sourceIndex] = (
                         f"Title: {source['title']}\n" +
                         f"Authors: {source['authors']}\n" +
                         f"URL: {source['url']}\n" +
                         f"Type: {source['type']}\n"
                         )
        # Loop over all associated test series to add nodes and edges
        for testType in source['tests']:
            if testType in ['vapor pressures','solubility limits','heat capacities']:
                for series in source['tests'][testType]:
                    seriesIndex += 1
                    nExpInSource += 1
                    currentSeries = source['tests'][testType][series]
                    # Calculate test series node positions
                    y[seriesIndex] = level2
                    x[seriesIndex] = 0
                    x[sourceIndex] += x[seriesIndex]
                    size[seriesIndex] = 20
                    names[seriesIndex] = series
                    types[seriesIndex] = testType
                    details[seriesIndex] = f"Status: {currentSeries['series status']}"
                    fullDetails[seriesIndex] = (
                                     f"Name: {series}\n" +
                                     f"Type: {testType}\n" +
                                     f"Status: {currentSeries['series status']}\n" +
                                     f"Database: {currentSeries['database']}\n" +
                                     f"Composition: {''.join([f'{key}: {currentSeries[compString][key]} ' for key in currentSeries[compString].keys()])}"
                                     )
                    # Create connections from sources to test series
                    edgeStarts.append(sourceIndex)
                    edgeEnds.append(seriesIndex)
                    # Get series statuses to set colors
                    seriesStatus = currentSeries['series status']
                    if seriesStatus == 'pass':
                        colorCode[seriesIndex] = RGB(0,255,0)
                    elif seriesStatus == 'fail':
                        colorCode[seriesIndex] = RGB(255,0,0)
                    elif seriesStatus == 'partial':
                        colorCode[seriesIndex] = RGB(255,255,0)
                    # Create sample points
                    nSampleInSeries = 0
                    for sample in currentSeries['samples'].keys():
                        sampleIndex += 1
                        nSampleInSeries += 1
                        currentSample = currentSeries['samples'][sample]
                        y[sampleIndex] = level3
                        x[sampleIndex] = samplePos
                        names[sampleIndex] = sample
                        types[sampleIndex] = testType
                        details[sampleIndex] = f"Status: {currentSample['status']}"
                        fullDetails[sampleIndex] = (
                                         f"Name: {series} sample {sample}\n" +
                                         f"Type: {testType}\n" +
                                         f"Status: {currentSample['status']}\n" +
                                         f"Database: {currentSeries['database']}\n" +
                                         f"Composition: {''.join([f'{key}: {currentSeries[compString][key]} ' for key in currentSeries[compString].keys()])}\n" +
                                         f"Temperature: {currentSample['temperature']}\n"
                                         )
                        if testType == 'vapor pressures':
                            fullDetails[sampleIndex] += f"Vapor pressures:\n{nlString.join([f'{key}: {currentSample[ppString][key]} atm' for key in currentSample[ppString].keys()])}\n"
                        elif testType == 'solubility limits':
                            fullDetails[sampleIndex] += f"Solubility limits:\n{nlString.join([f'{key}: {currentSample[frString][key]}' for key in currentSample[frString].keys()])}\n"
                        elif testType == 'heat capacities':
                            fullDetails[sampleIndex] += f"Heat capacity: {currentSample['heat capacity']} J/mol.K"
                        samplePos += sampleSpace
                        x[seriesIndex] += x[sampleIndex]
                        # Create connections from series to test samples
                        edgeStarts.append(seriesIndex)
                        edgeEnds.append(sampleIndex)
                        if currentSample['status'] == 'pass':
                            colorCode[sampleIndex] = RGB(0,255,0)
                        elif currentSample['status'] == 'fail':
                            colorCode[sampleIndex] = RGB(255,0,0)
                    x[seriesIndex] = x[seriesIndex] / nSampleInSeries
                    x[sourceIndex] += x[seriesIndex]
                    samplePos += sampleSpace
            elif testType == 'phase transitions':
                for series in source['tests'][testType]:
                    seriesIndex += 1
                    nExpInSource += 1
                    currentSeries = source['tests'][testType][series]
                    # Calculate test series node positions
                    y[seriesIndex] = level2
                    x[seriesIndex] = samplePos
                    x[sourceIndex] += x[seriesIndex]
                    samplePos += sampleSpace
                    names[seriesIndex] = series
                    types[seriesIndex] = testType
                    details[seriesIndex] = f"Status: {currentSeries['series status']}"
                    fullDetails[seriesIndex] = (
                                     f"Name: {series}\n" +
                                     f"Type: {testType}\n" +
                                     f"Status: {currentSeries['series status']}\n" +
                                     f"Database: {currentSeries['database']}\n" +
                                     f"Composition: {''.join([f'{key}: {currentSeries[compString][key]} ' for key in currentSeries[compString].keys()])}\n" +
                                     f"Target temperature: {currentSeries['target temperature']}K\n" +
                                     f"Low temperature phases: {' '.join(currentSeries['phases low'])}\n" +
                                     f"High temperature phases: {' '.join(currentSeries['phases high'])}\n" +
                                     f"Error tolerance: {currentSeries['error']}K\n"
                                     )
                    # Create connections from sources to test series
                    edgeStarts.append(sourceIndex)
                    edgeEnds.append(seriesIndex)
                    # Get series statuses to set colors
                    seriesStatus = currentSeries['series status']
                    if seriesStatus == 'pass':
                        colorCode[seriesIndex] = RGB(0,255,0)
                    elif seriesStatus == 'fail':
                        colorCode[seriesIndex] = RGB(255,0,0)
        x[sourceIndex] = x[sourceIndex] / nExpInSource

    graph = GraphRenderer()

    # Create nodes
    graph.node_renderer.glyph = Circle(size="size", fill_color="fill_color")
    graph.node_renderer.data_source.data = dict(index = node_indices,
                                                fill_color = colorCode,
                                                details = details,
                                                fullDetails = fullDetails,
                                                names = names,
                                                types = types,
                                                size = size)

    # Set hover and selection rules
    graph.node_renderer.selection_glyph = Circle(size=25, fill_color="fill_color")
    graph.node_renderer.hover_glyph = Circle(size=30, fill_color="fill_color")
    node_hover_tool = HoverTool(tooltips=[("Name", "@names"),("Type", "@types"),("Details", "@details")])
    graph.edge_renderer.selection_glyph = MultiLine(line_width=5)
    graph.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

    graph.selection_policy = NodesAndLinkedEdges()
    graph.inspection_policy = NodesAndLinkedEdges()

    plot = figure(title="MSTDB Tests", x_range=(-1,1), y_range=(-1,1),
                  tools=["xpan","xzoom_in","xzoom_out"], toolbar_location="right", plot_width=1800, plot_height=600)
    plot.axis.visible = False
    plot.xgrid.visible = False
    plot.ygrid.visible = False

    plot.add_tools(node_hover_tool, TapTool())

    graph_layout = dict(zip(node_indices, zip(x, y)))
    graph.edge_renderer.data_source.data = dict(start=edgeStarts,end=edgeEnds)
    graph.layout_provider = StaticLayoutProvider(graph_layout=graph_layout)

    plot.renderers.append(graph)

    # Create callback event to populate info panel for clicked node
    def nodeCallback(event):
        try:
            indexActive = graph.node_renderer.data_source.selected.indices[0]
            layout.children[2]=Div(text="<pre>" + fullDetails[indexActive] + "</pre>")
        except IndexError:
            pass

    plot.on_event(Tap, nodeCallback)

    return plot

# Filter buttons
sampleStatusOptions = ["pass", "fail", "incomplete"]
seriesStatusOptions = ["pass", "fail", "partial", "incomplete"]
seriesTypeOptions   = ["phase transitions", "solubility limits", "vapor pressures", "heat capacities"]
databaseOptions     = ["fluoride", "chloride"]
elementOptions      = ["Pu","U","Th","Nd","Ce","La","Cs","Zr","Rb","Ni","Fe","Cr","Ca","K","Al","Mg","Na","Be","Li","Cl","F"]
elementOptions.reverse()

def sampleStatusCallback(active):
    data.sampleStatusFilter = []
    for option in active:
        data.sampleStatusFilter.append(sampleStatusOptions[option])
    data.filter()
    plot = makeNetwork()
    layout.children[1] = plot

def seriesStatusCallback(active):
    data.seriesStatusFilter = []
    for option in active:
        data.seriesStatusFilter.append(seriesStatusOptions[option])
    data.filter()
    plot = makeNetwork()
    layout.children[1] = plot

def seriesTypeCallback(active):
    data.seriesTypeFilter = []
    for option in active:
        data.seriesTypeFilter.append(seriesTypeOptions[option])
    data.filter()
    plot = makeNetwork()
    layout.children[1] = plot

def seriesDatabaseCallback(active):
    data.seriesDatabaseFilter = []
    for option in active:
        data.seriesDatabaseFilter.append(databaseOptions[option])
    data.filter()
    plot = makeNetwork()
    layout.children[1] = plot

def seriesElementsCallback(active):
    data.seriesElementsFilter = []
    for option in active:
        data.seriesElementsFilter.append(elementOptions[option])
    data.filter()
    plot = makeNetwork()
    layout.children[1] = plot

sampleStatusButtonGroup = CheckboxButtonGroup(labels=sampleStatusOptions, active=[], width = 10)
sampleStatusButtonGroup.on_click(sampleStatusCallback)

seriesStatusButtonGroup = CheckboxButtonGroup(labels=seriesStatusOptions, active=[], width = 10)
seriesStatusButtonGroup.on_click(seriesStatusCallback)

seriesTypeButtonGroup = CheckboxButtonGroup(labels=seriesTypeOptions, active=[], width = 10)
seriesTypeButtonGroup.on_click(seriesTypeCallback)

seriesDatabaseButtonGroup = CheckboxButtonGroup(labels=databaseOptions, active=[], width = 10)
seriesDatabaseButtonGroup.on_click(seriesDatabaseCallback)

seriesElementsButtonGroup = CheckboxButtonGroup(labels=elementOptions, active=[], width = 10)
seriesElementsButtonGroup.on_click(seriesElementsCallback)

buttonRow = Column(
                Row(
                    Column(Div(text='Sample Status', width = 200),sampleStatusButtonGroup),
                    Column(Div(text='Series Status', width = 250),seriesStatusButtonGroup),
                    Column(Div(text='Series Type', width = 450),seriesTypeButtonGroup),
                    Column(Div(text='Database'),seriesDatabaseButtonGroup)
                ),
                Column(Div(text='Elements'),seriesElementsButtonGroup)
            )

filename = 'verificationData-tested.json'
data = parseTests.jsonTestData(filename)

div=Div(text='')

plot = makeNetwork()

layout = Column(buttonRow,plot,div)

curdoc().add_root(layout)

# output_file("interactive_graphs.html")
# show(plot)
