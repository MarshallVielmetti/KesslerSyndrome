using SatelliteToolbox

function propogate(TLE)
    orbp = init_orbit_propagator(Val(:sgp4), TLE)
    r, v = propagate!(orbp, (0:3:24)*60*60)
    return r, v
end

function createTLEItems()

end


function runSimulation() 
    init_propogator()
    TLE = createTLEItems()
    rs, vs = propogate.(TLE)
end