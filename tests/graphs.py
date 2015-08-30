import plotly.plotly as py
from plotly.graph_objs import *

def box():
    trace0 = Box(
        y = [
        11.743771306564682,
        0.008477790012827552,
        1.3601017526409156e-05,
        0.06316186780728096,
        4.509241323978885e-05,
        0.0076751327330669165,
        0.1407651271022885
        ]
    )
    trace1 = Box(
        y = [
        11.743771306564682,
        0.008477790012827552,
        1.3601017526409156e-05,
        0.04688244882832312,
        4.509241323978885e-05,
        23.02651565202709,
        3.939861080619553e-05,
        0.0076751327330669165,
        2.859234616963753
        ]
    )
    trace2 = Box(
        y = [
        11.743771306564682, 
        0.008477790012827552, 
        0.9260989586841414, 
        2.5545608699467885e-05, 
        1.3601017526409156e-05, 
        0.04688244882832312, 
        4.509241323978885e-05, 
        23.02651565202709, 
        1336.1235651108411, 
        24.670583045393737, 
        3.939861080619553e-05, 
        0.0076751327330669165, 
        2.859234616963753
        ]
    )
    data = Data([trace0, trace1, trace2])
    plot_url = py.plot(data, filename='distance-boxes')

def scatter():
    trace0 = Scatter(
        x=[
        20.786516853932586,
        21.910112359550563,
        42.3728813559322,
        51.724137931034484,
        58.620689655172406,
        19.54022988505747,
        63.2183908045977,
        73.83720930232558,
        39.53488372093023,
        27.325581395348834
        ],
        y=[
        100.0,
        100.0,
        85.71428571428571,
        77.77777777777779,
        60,
        50,
        69.23076923076923,
        57.14285714285714,
        85.71428571428571,
        80
        ],
        mode='markers+text',
        text=[
        '0', # none
        '1', # no email
        '2', # dept removed
        '3', # First line removed
        '4', # All but first line, department removed
        '5', # Last 2 lines
        '6', # Second and third lines only
        '7', # Second and last lines 
        '8', # Last 3 lines
        '9'  # All but first 2 lines
        ],
        textposition='top center',
        textfont=Font(
            family='Arial',
        ),
        marker=Marker(
            size=12,
        )
    )

    data = Data([trace0])
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

if __name__ == "__main__":
    box()