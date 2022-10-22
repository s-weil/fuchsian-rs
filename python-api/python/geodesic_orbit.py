import python_api


horocyclic = python_api.moebius_matrix(1.0, 1.0, 0.0, 1.0)
print(horocyclic)

elliptic= python_api.moebius_matrix(0.0, -1.0, 1.0, 0.0)
print(elliptic)

modular_group = [horocyclic, elliptic]

geodesic_end_points = (-1.0, 1.0) 
# geodesic_orbit = python_api.geodesic_orbit(modular_group, geodesic_end_points, 10_000, 20, "random")
geodesic_orbit = python_api.geodesic_orbit(modular_group, geodesic_end_points, 1_000, 20, "sequential")
print("generated geodesic orbit")


###### PLOTTING ######

import plotly.graph_objects as go

fig = go.Figure()

for geodesic in geodesic_orbit:
    x = [ z[0] for z in geodesic ]
    y = [ z[1] for z in geodesic ]
    fig.add_trace(go.Line(x=x, y=y, mode='lines'))

fig.update_layout(
    title="gedeosic orbits", xaxis_title="Re", yaxis_title="Im"
)
    
fig.show()