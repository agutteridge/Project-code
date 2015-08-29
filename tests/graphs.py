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

if __name__ == "__main__":
	trace0 = Scatter(
	    x=[
	    20.786516853932586,
	    21.910112359550563,
	    42.3728813559322,
	    77.96610169491525
	    ],
	    y=[
	    100.0,
	    100.0,
	    85.71428571428571,
	    23.52941176470588
	    ],
	    mode='markers+text',
	    name='At least 1 line',
	    text=[
	    'n/a, email removed',
	    ' ',
	    'Department removed',
	    'First'
	    ],
	    textposition='top center',
	    textfont=Font(
	        family='Arial',
	    ),
	    marker=Marker(
	        size=12,
	    )
	)
	trace1 = Scatter(
	    x=[
	    52.87356321839081,
	    59.32203389830508,
	    67.81609195402298,
	    19.54022988505747
	    ],
	    y=[
	    45.45454545454545,
	    60,
	    40,
	    50.0
	    ],
	    mode='markers+text',
	    name='At least 2 lines',
	    text=[
	    'Second',
	    'First two',
	    'First and last',
	    'Last two',
	    ],
	    textposition='top center',
	    textfont=Font(
	        family='Arial',
	    ),
	    marker=Marker(
	        size=12,
	    )
	)
	trace2 = Scatter(
	    x=[
	    44.44444444444444,
	    33.33333333333333,
	    39.76608187134503,
	    51.461988304093566
	    ],
	    y=[
	    75.0,
	    71.42857142857143,
	    85.71428571428571,
	    62.5
	    ],
	    mode='markers+text',
	    name='At least 3 lines',
	    text=[
	    'Second and last',
	    'Third',
	    'Last three',
	    'First three'
	    ],
	    textposition='top center',
	    textfont=Font(
	        family='Arial',
	    ),
	    marker=Marker(
	        size=12,
	    )
	)

	trace3 = Scatter(
	    x=[
	    	41.221374045801525,
	    	31.297709923664126
	    ],
	    y=[
	    	100,
	    	100
	    ],
	    mode='markers+text',
	    name='At least 4 lines',
	    text=[
	    'First 4',
	    'Last 4'
	    ],
	    textposition='top center',
	    textfont=Font(
	        family='Arial',
	    ),
	    marker=Marker(
	        size=12,
	    )
	)

	data = Data([trace0, trace1, trace2, trace3])
	layout = Layout(
	    xaxis=XAxis(
	        autorange=True,
	    ),
	    yaxis=YAxis(
	        autorange=True,
	    ),
	    legend=Legend(
	        y=0.5,
	        yref='paper',
	        font=Font(
	            size=18,
	        )
	    ),
	)

	fig = Figure(data=data, layout=layout)
	plot_url = py.plot(fig, filename='geocode_performance_scatter')