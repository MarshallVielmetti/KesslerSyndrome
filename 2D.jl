using Plots
using Random

#For each 'object' I want to define an ( x,y ), (x,y) velocity, and 'mass'

mutable struct Obj
    xpos::Float32
    ypos::Float32
    xvel::Float32
    yvel::Float32
    mass::Float32
    radius::Float32
end

#Obj() = Obj()

rn = rand(Float64)
print(rn)