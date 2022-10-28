import python_api
import plotly.graph_objects as go

###########################
###### MODULAR GROUP ######

horocyclic = python_api.moebius_matrix(1.0, 1.0, 0.0, 1.0)
print(horocyclic)

elliptic= python_api.moebius_matrix(0.0, -1.0, 1.0, 0.0)
print(elliptic)

# modular_group = [horocyclic, elliptic]
modular_group = [horocyclic, elliptic, elliptic, horocyclic]

# Euclidean height
horocycle_height = 2.0
horocyclic_orbit = python_api.horocyclic_orbit(modular_group, horocycle_height, 100, 30, "sequential")
print("generated horocyclic orbit")


###### PLOTTING ######

fig = go.Figure()

for horocycle in horocyclic_orbit:
    x = [ z[0] for z in horocycle ]
    y = [ z[1] for z in horocycle ]
    fig.add_trace(go.Line(x=x, y=y, mode='lines'))

fig.update_layout(
    title="horocycle orbits", xaxis_title="Re", yaxis_title="Im"
)
    
fig.show()


###########################
###### HYPERBOLIC #########

# hyperbolic has 2 fixed points on the boundary, 0 and infty
hyperbolic = python_api.moebius_matrix(5.0, 0.0, 0.0, 0.2)
single_group = [ hyperbolic ]

# Euclidean height
horocycle_height = 1.0
horocyclic_orbit = python_api.horocyclic_orbit(single_group, horocycle_height, 100, 30, "sequential")
print("generated horocyclic orbit")

###### PLOTTING ######

fig = go.Figure()

for horocycle in horocyclic_orbit:
    x = [ z[0] for z in horocycle ]
    y = [ z[1] for z in horocycle ]
    fig.add_trace(go.Line(x=x, y=y, mode='lines'))

fig.update_layout(
    title="horocycle orbits", xaxis_title="Re", yaxis_title="Im"
)
    
fig.show()
