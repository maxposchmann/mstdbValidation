from bokeh.io import show, output_file
from bokeh.models import (Circle, MultiLine,
                          GraphRenderer, StaticLayoutProvider, NodesAndLinkedEdges,
                          Div, Column,
                          HoverTool, TapTool, BoxSelectTool)
from bokeh.layouts import column
from bokeh.plotting import figure, curdoc
from bokeh.colors import RGB
from bokeh.palettes import Spectral4
from bokeh.events import Tap
import parseTests

filename = 'verificationData-results.json'
data = parseTests.jsonTestData(filename)

totalNodes = len(data.experimentalReferences) + len(data.testSeries)
node_indices = list(range(totalNodes))

seriesNames = [list(series.keys())[0] for series in data.testSeries]

x = []
y = []
names = []
types = []
details = []
fullDetails = []
size = []
colorCode = [RGB(100,100,100) for i in node_indices]
sizes = [50 for i in node_indices]

level1 = 0.5
level2 = -0.5
leftLim = -0.75
rightLim = 0.75

edgeStarts = []
edgeEnds   = []
i = -1
j = -1
for source in data.experimentalReferences:
    i += 1
    # Calculate node reference positions
    y.append(level1)
    x.append(leftLim + (rightLim - leftLim) * (i+0.5) / len(data.experimentalReferences))
    colorCode[i] = RGB(0,0,255)
    size.append(20)
    names.append(source['title'])
    types.append(source['type'])
    details.append(f"Authors: {source['authors']}")
    fullDetails.append("<pre>" +
                       f"Title: {source['title']}\n" +
                       f"Authors: {source['authors']}\n" +
                       f"URL: {source['url']}\n" +
                       f"Type: {source['type']}\n" +
                       "</pre>")
    # Loop over all associated test series to add nodes and edges
    for testType in source['tests']:
        if testType == 'vapor pressures':
            for series in source['tests'][testType]:
                j += 1
                currentSeries = source['tests'][testType][series]
                # Calculate test series node positions
                y.append(level2)
                x.append(leftLim + (rightLim - leftLim) * (j+0.5) / len(data.testSeries))
                size.append(20)
                names.append(series)
                types.append(currentSeries['series type'])
                details.append(f"Status: {currentSeries['series status']}")
                tempTest = 'tests'
                tempComp = 'composition'
                fullDetails.append("<pre>" +
                                   f"Name: {series}\n" +
                                   f"Type: {currentSeries['series type']}\n" +
                                   f"Status: {currentSeries['series status']}\n" +
                                   f"Database: {currentSeries['database']}\n" +
                                   f"Composition: {''.join([f'{key}: {currentSeries[tempComp][key]} ' for key in currentSeries['composition'].keys()])}" +
                                   "</pre>")
                # Create connections from sources to test series
                edgeStarts.append(i)
                edgeEnds.append(j + len(data.experimentalReferences))
                # Get series statuses to set colors
                seriesStatus = currentSeries['series status']
                if seriesStatus == 'pass':
                    colorCode[j + len(data.experimentalReferences)] = RGB(0,255,0)
                elif seriesStatus == 'fail':
                    colorCode[j + len(data.experimentalReferences)] = RGB(255,0,0)
                elif seriesStatus == 'partial':
                    colorCode[j + len(data.experimentalReferences)] = RGB(255,255,0)

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

plot = figure(title="MSTDB Tests", x_range=(-1.1,1.1), y_range=(-1.1,1.1),
              tools="", toolbar_location=None)
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
        layout.children[1]=Div(text=fullDetails[indexActive])
    except IndexError:
        pass

plot.on_event(Tap, nodeCallback)

div=Div(text='')

layout=Column(plot,div)
curdoc().add_root(layout)

# output_file("interactive_graphs.html")
# show(plot)
