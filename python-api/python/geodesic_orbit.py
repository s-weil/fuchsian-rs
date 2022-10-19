import python_api


horocyclic = python_api.moebius_matrix(1.0, 1.0, 0.0, 1.0)
print(horocyclic)

elliptic= python_api.moebius_matrix(0.0, -1.0, 1.0, 0.0)
print(elliptic)

modular_group = [horocyclic, elliptic]

geodesic_end_points = (-1.0, 1.0) # two fixed points
geodesic_orbit = python_api.geodesic_orbit(modular_group, geodesic_end_points, 100, 50, "random")
print("generated geodesic orbit")


###### PLOTTING ######

import plotly.graph_objects as px

TODO
