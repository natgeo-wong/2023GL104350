using DrWatson
@quickactivate "2023GL104350"

include(srcdir("extractprecip.jl"))

extractprecip("DGW","S1284km300V64","damping001")
