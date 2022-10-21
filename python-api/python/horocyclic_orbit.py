import python_api


horocyclic = python_api.moebius_matrix(1.0, 1.0, 0.0, 1.0)
print(horocyclic)

elliptic= python_api.moebius_matrix(0.0, -1.0, 1.0, 0.0)
print(elliptic)

modular_group = [horocyclic, elliptic]

horocycle_height = 2.0
horocyclic_orbit = python_api.horocyclic_orbit(modular_group, horocycle_height, 10, 30, "random")
print("generated geodesic orbit")


###### PLOTTING ######

import plotly.graph_objects as go

fig = go.Figure()

for horocycle in horocyclic_orbit:
    x = [ z[0] for z in horocycle ]
    y = [ z[1] for z in horocycle ]
    fig.add_trace(go.Line(x=x, y=y, mode='lines'))

fig.update_layout(
    title="horocycle orbits", xaxis_title="Re", yaxis_title="Im"
)
    
fig.show()