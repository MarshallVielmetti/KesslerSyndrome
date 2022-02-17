using SatelliteToolbox
import Random
using Plots
using StatsBase
using LaTeXStrings
gr()

global eop = get_iers_eop()

# The DCM (Direction Cosine Matrix) that rotates TEME into alignment with ITRF
global D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), DatetoJD(2022,1,1,0,0,0), eop)

function runAnalysis()
    tles = read_tle("TLEData.txt")

    # Random sample of 3000 TLEs
    ids = floor.(rand((1:length(tles)), 3000))

    rand_tles = [tles[x] for (x) in ids]

    altitudes = getTleAltitude.(tles)
    filter!((x) -> x.<2500 && x.>100, altitudes)

    return altitudes
end

function getTleAltitude(tle)
    orbp = init_orbit_propagator(Val(:sgp4), tle)
    try
        r_teme, v_teme = propagate_to_epoch!(orbp, DatetoJD(2022, 1, 12, 0, 0, 0))
        r_itrf = D_ITRF_TEME * r_teme
        lat,lon,h = ecef_to_geodetic(r_itrf)
    
        return h / 1000;
    catch
        return 0
    end
end

function doVisualization(altitudes)
    histogram(altitudes, bins=range(0,step=50,stop=2500), label="Tracked Objects", title="Proportion of Space Debris by Orbital Altitude", color="azure2", normalize=true)
    xlabel!("Altitude (km)")
    ylabel!("Proportion of Debris")
end

function doSpatialDensityCalculations(altitudes)
    bins = 150:25:2500

    h = fit(Histogram, altitudes, bins)

    bins = collect(bins)
    popat!(bins, 1)

    densities = getSphereDensity.(bins, h.weights)

    plot(bins, densities, color="black", label="Spatial Density", yscale = :log10, ylims = (10E-10, 10E-4))
    title!("Spatial Density by Altitude")
    xlabel!(L"Altitude ($km$)")
    ylabel!(L"Spatial Density ($No. / Km^3$)")
end

function getSphereDensity(radius, count)
    outerVolume = 4/3 * π * radius^3
    innerVolume = 4/3 * pi * (radius - 5)^3

    sliceVolume = outerVolume - innerVolume

    return count / sliceVolume
end

function doInclinationCalculation() 

end



altitudes = runAnalysis()
doVisualization(altitudes)
doSpatialDensityCalculations(altitudes)


png("Spatial Density")