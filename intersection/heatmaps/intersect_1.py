import plotly as py
import plotly.graph_objs as go

trace = go.Heatmap(
	z=[[0.8, 0.61, 0.608, 0.365, 0.486, 0.411, 0.441, 0.431, 0.514], [0.61, 0.8, 0.611, 0.439, 0.422, 0.405, 0.55, 0.492, 0.525], [0.608, 0.611, 0.8, 0.431, 0.592, 0.498, 0.548, 0.505, 0.548], [0.365, 0.439, 0.431, 0.8, 0.417, 0.468, 0.536, 0.506, 0.426], [0.486, 0.422, 0.592, 0.417, 0.8, 0.515, 0.456, 0.5, 0.552], [0.411, 0.405, 0.498, 0.468, 0.515, 0.8, 0.423, 0.426, 0.489], [0.441, 0.55, 0.548, 0.536, 0.456, 0.423, 0.8, 0.678, 0.5], [0.431, 0.492, 0.505, 0.506, 0.5, 0.426, 0.678, 0.8, 0.513], [0.514, 0.525, 0.548, 0.426, 0.552, 0.489, 0.5, 0.513, 0.8]],
	x=["D206_m2", "D206_m9", "D206_mR", "D410_m1", 
	"D410_m5", "D410_mR", "D430_m7", "D430_m8", "D430_mR"],
	y=["D206_m2", "D206_m9", "D206_mR", "D410_m1", 
	"D410_m5", "D410_mR", "D430_m7", "D430_m8", "D430_mR"])
data=[trace]
py.offline.plot(data, filename='heatmap_1')