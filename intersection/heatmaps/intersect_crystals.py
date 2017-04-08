import plotly as py
import plotly.graph_objs as go

trace = go.Heatmap(
	z=[[0.6, 0.498, 0.548], [0.498, 0.6, 0.489], [0.548, 0.489, 0.6]],
	x=["D206_mR", "D410_mR", "D430_mR"],
	y=["D206_mR", "D410_mR", "D430_mR"])
data=[trace]
py.offline.plot(data, filename='heatmap_c')