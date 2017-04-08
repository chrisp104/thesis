import plotly as py
import plotly.graph_objs as go

trace = go.Heatmap(
	z=[[33, 14, 13, 10, 11, 11, 12, 12, 13], [14, 33, 15, 13, 11, 11, 16, 13, 12], [13, 15, 33, 10, 13, 13, 12, 10, 12], [10, 13, 10, 33, 10, 13, 15, 13, 10], [11, 11, 13, 10, 33, 12, 11, 12, 14], [11, 11, 13, 13, 12, 33, 11, 12, 13], [12, 16, 12, 15, 11, 11, 33, 16, 12], [12, 13, 10, 13, 12, 12, 16, 33, 13], [13, 12, 12, 10, 14, 13, 12, 13, 33]],
	x=["D206_m2", "D206_m9", "D206_mR", "D410_m1", 
	"D410_m5", "D410_mR", "D430_m7", "D430_m8", "D430_mR"],
	y=["D206_m2", "D206_m9", "D206_mR", "D410_m1", 
	"D410_m5", "D410_mR", "D430_m7", "D430_m8", "D430_mR"])
data=[trace]
py.offline.plot(data, filename='basic-heatmap')