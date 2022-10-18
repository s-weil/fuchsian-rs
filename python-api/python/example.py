import python_api


horocyclic = python_api.moebius_matrix(1.0, 1.0, 0.0, 1.0);
print(horocyclic)

elliptic= python_api.moebius_matrix(0.0, -1.0, 1.0, 0.0);
print(elliptic)

modular_group = [horocyclic, elliptic]

base_point = (0.5, 1.5) # Note that i is a singularity of the modular group
orbit = python_api.orbit(modular_group, base_point, 100)
print(orbit)