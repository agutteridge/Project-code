import plotly.plotly as py
from plotly.graph_objs import *

def plot(x):
	data = Data([
		Histogram(
			x=x
		)
	])

	layout = Layout(
	    yaxis=YAxis(
	        type='log',
	        autorange=True
	    )
	)
	plot_url = py.plot(data, filename='address-lines-histogram')