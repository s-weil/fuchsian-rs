import python_api as fuchsian

horocyclic = fuchsian.moebius_matrix(1.0, 1.0, 0.0, 1.0)
print(horocyclic)

elliptic= fuchsian.moebius_matrix(0.0, -1.0, 1.0, 0.0)
print(elliptic)

modular_group = [horocyclic, elliptic]

base_point = (1.0, 0.0) # Note that i = (0, 1) is a singularity of the modular group
# orbit = fuchsian.orbit(modular_group, base_point, 100_000, "random")
orbit = fuchsian.orbit(modular_group, base_point, 100_000, "sequential")
print("generated orbit")


###### PLOTTING ######

import plotly.graph_objects as px

x = [ z[0] for z in orbit ]
y = [ z[1] for z in orbit ]

plot = px.Figure(data=[px.Scatter(
	x=x,
	y=y,
	mode='markers')
])
plot.show()
