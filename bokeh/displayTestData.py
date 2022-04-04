from bokeh.io import show, output_file
from bokeh.models import (Circle, MultiLine,
                          GraphRenderer, StaticLayoutProvider, NodesAndLinkedEdges,
                          Div, Column, Row,
                          HoverTool, TapTool, BoxSelectTool,
                          CheckboxButtonGroup, CustomJS)
from bokeh.layouts import column
from bokeh.plotting import figure, curdoc
from bokeh.colors import RGB
from bokeh.palettes import Spectral4
from bokeh.events import Tap
import parseTests
import testRunner

def makeNetwork():
    nExpRef = data.nSources
    nTestSeries = data.nSeries
    totalNodes = nExpRef + nTestSeries
    node_indices = list(range(totalNodes))

    x = [0 for i in node_indices]
    y = [0 for i in node_indices]
    names = ['' for i in node_indices]
    types = ['' for i in node_indices]
    details = ['' for i in node_indices]
    fullDetails = ['' for i in node_indices]
    size = [0 for i in node_indices]
    colorCode = [RGB(100,100,100) for i in node_indices]
    sizes = [50 for i in node_indices]

    level1 = 0.5
    level2 = -0.5
    leftLim = -0.95
    rightLim = 0.95

    edgeStarts = []
    edgeEnds   = []
    i = -1
    j = nExpRef - 1
    for sourceName in data.data['sources']:
        source = data.data['sources'][sourceName]
        i += 1
        # Calculate node reference positions
        y[i] = level1
        x[i] = 0
        nExpInSource = 0
        colorCode[i] = RGB(0,0,255)
        size[i] = 20
        names[i] = source['title']
        types[i] = source['type']
        details[i] = f"Authors: {source['authors']}"
        fullDetails[i] = (
                         f"Title: {source['title']}\n" +
                         f"Authors: {source['authors']}\n" +
                         f"URL: {source['url']}\n" +
                         f"Type: {source['type']}\n"
                         )
        # Loop over all associated test series to add nodes and edges
        for testType in source['tests']:
            if testType in ['vapor pressures','solubility limits','heat capacities']:
                for series in source['tests'][testType]:
                    j += 1
                    nExpInSource += 1
                    currentSeries = source['tests'][testType][series]
                    # Calculate test series node positions
                    y[j] = level2
                    x[j] = leftLim + (rightLim - leftLim) * (j-nExpRef+0.5) / nTestSeries
                    x[i] += x[j]
                    size[j] = 20
                    names[j] = series
                    types[j] = testType
                    details[j] = f"Status: {currentSeries['series status']}"
                    tempComp = 'composition'
                    fullDetails[j] = (
                                     f"Name: {series}\n" +
                                     f"Type: {testType}\n" +
                                     f"Status: {currentSeries['series status']}\n" +
                                     f"Database: {currentSeries['database']}\n" +
                                     f"Composition: {''.join([f'{key}: {currentSeries[tempComp][key]} ' for key in currentSeries['composition'].keys()])}"
                                     )
                    # Create connections from sources to test series
                    edgeStarts.append(i)
                    edgeEnds.append(j)
                    # Get series statuses to set colors
                    seriesStatus = currentSeries['series status']
                    if seriesStatus == 'pass':
                        colorCode[j] = RGB(0,255,0)
                    elif seriesStatus == 'fail':
                        colorCode[j] = RGB(255,0,0)
                    elif seriesStatus == 'partial':
                        colorCode[j] = RGB(255,255,0)
            elif testType == 'phase transitions':
                for series in source['tests'][testType]:
                    j += 1
                    nExpInSource += 1
                    currentSeries = source['tests'][testType][series]
                    # Calculate test series node positions
                    y[j] = level2
                    x[j] = leftLim + (rightLim - leftLim) * (j-nExpRef+0.5) / nTestSeries
                    x[i] += x[j]
                    size[j] = 20
                    names[j] = series
                    types[j] = testType
                    details[j] = f"Status: {currentSeries['series status']}"
                    tempComp = 'composition'
                    fullDetails[j] = (
                                     f"Name: {series}\n" +
                                     f"Type: {testType}\n" +
                                     f"Status: {currentSeries['series status']}\n" +
                                     f"Database: {currentSeries['database']}\n" +
                                     f"Composition: {''.join([f'{key}: {currentSeries[tempComp][key]} ' for key in currentSeries['composition'].keys()])}\n" +
                                     f"Target temperature: {currentSeries['target temperature']}K\n" +
                                     f"Low temperature phases: {' '.join(currentSeries['phases low'])}\n" +
                                     f"High temperature phases: {' '.join(currentSeries['phases high'])}\n" +
                                     f"Error tolerance: {currentSeries['error']}K\n"
                                     )
                    # Create connections from sources to test series
                    edgeStarts.append(i)
                    edgeEnds.append(j)
                    # Get series statuses to set colors
                    seriesStatus = currentSeries['series status']
                    if seriesStatus == 'pass':
                        colorCode[j] = RGB(0,255,0)
                    elif seriesStatus == 'fail':
                        colorCode[j] = RGB(255,0,0)
        x[i] = x[i] / nExpInSource

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
                  tools="", toolbar_location=None, plot_width=1800, plot_height=600)
    plot.axis.visible = False
    plot.xgrid.visible = False
    plot.ygrid.visible = False

    plot.add_tools(node_hover_tool, TapTool(), BoxSelectTool())

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

seriesElementsButtonGroup = CheckboxButtonGroup(labels=elementOptions, active=[], width = 10)
seriesElementsButtonGroup.on_click(seriesElementsCallback)

buttonRow = Column(
                Row(
                    Column(Div(text='Sample Status', width = 200),sampleStatusButtonGroup),
                    Column(Div(text='Series Status', width = 250),seriesStatusButtonGroup),
                    Column(Div(text='Series Type'),seriesTypeButtonGroup)
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
