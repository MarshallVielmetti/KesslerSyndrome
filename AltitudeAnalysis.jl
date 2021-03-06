using SatelliteToolbox
import Random
using Plots
using StatsBase
using LaTeXStrings
gr()

global eop = get_iers_eop()

# The DCM (Direction Cosine Matrix) that rotates TEME into alignment with ITRF
global D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), DatetoJD(2022, 1, 1, 0, 0, 0), eop)

function runAnalysis()
    tles = read_tle("TLEData.txt")

    altitudes = getTleAltitude.(tles)
    filter!((x) -> x .< 2500 && x .> 100, altitudes)

    return altitudes
end

function getTleAltitude(tle)
    orbp = init_orbit_propagator(Val(:sgp4), tle)
    try
        r_teme, v_teme = propagate_to_epoch!(orbp, DatetoJD(2022, 1, 12, 0, 0, 0))
        r_itrf = D_ITRF_TEME * r_teme
        lat, lon, h = ecef_to_geodetic(r_itrf)

        return h / 1000
    catch
        return 0
    end
end

# Returns tupple altitude, inclination
function getTleAltitudeInclination(tle)
    orbp = init_orbit_propagator(Val(:sgp4), tle)
    try
        r_teme, v_teme = propagate_to_epoch!(orbp, DatetoJD(2022, 1, 12, 0, 0, 0))
        r_itrf = D_ITRF_TEME * r_teme
        lat, lon, h = ecef_to_geodetic(r_itrf)

        return [h / 1000; tle.i]
    catch
        return [0; 0]
    end
end

function doVisualization(altitudes)
    histogram(altitudes, bins = range(0, step = 50, stop = 2500), label = "Tracked Objects", title = "Proportion of Space Debris by Orbital Altitude", color = "azure2", normalize = true)
    xlabel!("Altitude (km)")
    ylabel!("Proportion of Debris")
end

function doSpatialDensityCalculations(altitudes)
    bins = 150:25:2500

    h = fit(Histogram, altitudes, bins)

    bins = collect(bins)
    popat!(bins, 1)

    densities = getSphereDensity.(bins, h.weights)

    plot(bins, densities, color = "black", label = "Spatial Density", yscale = :log10, ylims = (10E-10, 10E-4))
    title!("Spatial Density by Altitude")
    xlabel!(L"Altitude ($km$)")
    ylabel!(L"Spatial Density ($No. / Km^3$)")
end

function getSphereDensity(radius, count)
    outerVolume = 4 / 3 * π * radius^3
    innerVolume = 4 / 3 * pi * (radius - 0.005)^3

    sliceVolume = outerVolume - innerVolume

    return count / sliceVolume
end

#Main driver function for inclination analysis
function doInclinationAnalysis()
    # Plot relative densities in altitude by inclination - 2D histogram
    # Using 1km and 1 degree bins
    tles = read_tle("TLEData.txt")

    # m = matrix [altitude, inclination]
    m = getTleAltitudeInclination.(tles)

    m = mapreduce(permutedims, vcat, m)

    m = m[(m[:, 1].<2500).&(m[:, 1].>100), :]

    altitudes = m[:, 1]
    inclinations = m[:, 2]

    # weights array altitudes by inc
    d = fit(Histogram, (altitudes, inclinations), ((100:10:2500, 0:1:130)), closed = :right)

    bindata = zeros(1, 3)

    for alt in 1:size(d.weights)[1]
        for incl in 1:size(d.weights)[2]
            if d.weights[alt, incl] > 0
                bindata = [bindata; [alt * 10 incl d.weights[alt, incl]]]
            end
        end
    end
    @show size(bindata)

    bindata = bindata[2:end, :]

    bindata[:, 3] = log.(bindata[:, 3])

    @show size(bindata)

    plot(bindata[:, 1], bindata[:, 2], seriestype = :scatter, markersize = bindata[:, 3], color = "lightblue", legend = false)
    title!("Altitude by Orbital Inclination, Count of Objects")
    xlabel!("Altitude (10km Bins)")
    ylabel!("Inclination, 1° Bins")
    png("Figures/PNG/Altitude by Orbital Inclination, Count of Objects")
    savefig("Figures/SVG/Altitude by Orbital Inclination, Count of Objects.svg")
    # histogram2d(m[:, 1], m[:, 2], bins = (100:10:2500, 0:1:130))
    # return altitudes
end

#Main driver function for inclination analysis for only deb objects
function doDebrisOnlyInclinationAnalysis()
    # Plot relative densities in altitude by inclination - 2D histogram
    # Using 1km and 1 degree bins
    tles = read_tle("TLEData.txt")

    # @show fieldnames(typeof(tles[1]))
    # @show tles[1].name
    # @show occursin("Deb", tles[1].name)

    filter!(e -> occursin("DEB", e.name), tles)

    # m = matrix [altitude, inclination]
    m = getTleAltitudeInclination.(tles)

    m = mapreduce(permutedims, vcat, m)

    m = m[(m[:, 1].<2500).&(m[:, 1].>100), :]

    altitudes = m[:, 1]
    inclinations = m[:, 2]

    # weights array altitudes by inc
    d = fit(Histogram, (altitudes, inclinations), ((100:10:2500, 0:1:130)), closed = :right)

    bindata = zeros(1, 3)

    for alt in 1:size(d.weights)[1]
        for incl in 1:size(d.weights)[2]
            if d.weights[alt, incl] > 0
                bindata = [bindata; [alt * 10 incl d.weights[alt, incl]]]
            end
        end
    end
    @show size(bindata)

    bindata = bindata[2:end, :]

    bindata[:, 3] = log.(bindata[:, 3])

    @show size(bindata)

    plot(bindata[:, 1], bindata[:, 2], seriestype = :scatter, markersize = bindata[:, 3], color = :darkred, legend = false)
    title!("Altitude by Orbital Inclination, Count of Debris")

    xlabel!("Altitude (10km Bins)")
    ylabel!("Inclination, 1° Bins")
    png("Figures/PNG/Debris - Orbit Inclination x Altitude Graph")
    savefig("Figures/SVG/Debris - Orbit Inclination x Altitude Graph.svg")
    # return altitudes
end

function doDebrisAltitudeBinAnalysis()
    tles = read_tle("TLEData.txt")

    filter!(e -> occursin("DEB", e.name), tles)

    altitudes = getTleAltitude.(tles)
    filter!((x) -> x .< 2500 && x .> 100, altitudes)

    histogram(altitudes, bins = range(0, step = 50, stop = 2500), label = "Tracked Objects", title = "Proportion of Space Debris by Orbital Altitude", color = "darkred", normalize = true)
    xlabel!("Altitude (km)")
    ylabel!("Proportion of Debris")

    png("Figures/PNG/Debris - Proportion by Altitude")
    savefig("Figures/SVG/Debris - Proportion by Altitude.svg")

end

function doDebrisSpacialDensityAnalysis()
    tles = read_tle("TLEData.txt")

    filter!(e -> occursin("DEB", e.name), tles)

    altitudes = getTleAltitude.(tles)
    filter!((x) -> x .< 2500 && x .> 100, altitudes)

    bins = 150:25:2500

    h = fit(Histogram, altitudes, bins)

    bins = collect(bins)
    popat!(bins, 1)

    densities = getSphereDensity.(bins, h.weights)

    plot(bins, densities, color = "black", label = "Spatial Density", yscale = :log10, ylims = (10E-10, 10E-4))
    title!("Spatial Density of Debris by Altitude")
    xlabel!(L"Altitude ($km$)")
    ylabel!(L"Spatial Density ($No. / Km^3$)")

    png("Figures/PNG/Debris - Spatial Density")
    savefig("Figures/SVG/Debris - Spatial Density.svg")
end

function getStatistics()
    tles = read_tle("TLEData.txt")

    @show length(tles)

    filter!(e -> occursin("DEB", e.name), tles)

    @show length(tles)
end

function doProportionDebrisByAltitude()
    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    filter!((x) -> x .< 2500 && x .> 100, DEB_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, NOD_altitudes)

    bins = (0:50:2500)

    DEB_binned = fit(Histogram, DEB_altitudes, bins, closed = :right)
    NOD_binned = fit(Histogram, NOD_altitudes, bins, closed = :right)

    bins_arr = collect(bins)
    popat!(bins_arr, 1)

    proportions = DEB_binned.weights ./ (DEB_binned.weights .+ NOD_binned.weights)
    plot(bins_arr, proportions, color = "black", legend = false, ylims = (0, 1))

    title!("Proportion of Objects Classified as Debris by Altitude")
    xlabel!(L"Altitude ($km$)")
    # xlims!(0, 1)
    ylabel!("Proportion of Objects")

    png("Figures/PNG/Comb - Proportion Debris Objects")
    savefig("Figures/SVG/Comb - Proportion Debris Objects.svg")

end

function doCountDebrisNonDebrisComparisonByAltitude()
    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    filter!((x) -> x .< 2500 && x .> 100, DEB_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, NOD_altitudes)

    bins = (0:50:2500)

    DEB_binned = fit(Histogram, DEB_altitudes, bins, closed = :right)
    NOD_binned = fit(Histogram, NOD_altitudes, bins, closed = :right)

    bins_arr = collect(bins)
    popat!(bins_arr, 1)

    DEB_ADJ = DEB_binned.weights
    NONF = NOD_binned.weights .* 0.35 #Non functional satellites
    NOD_ADJ = NOD_binned.weights .* (1 - 0.35) #Functional satellites

    VALS = [DEB_ADJ NONF NOD_ADJ]

    plot(bins_arr, VALS, color = ["darkred" "green" "lightblue"], legend = true, label = ["Debris" "Non-Functional Satellites" "Functional Satellites"])

    title!("Counts of Objects By Debris Status")
    xlabel!(L"Altitude ($km$)")
    # xlims!(0, 1)
    ylabel!("Count of Objects")

    png("Figures/PNG/Comb - Count of Objects By Debris Status")
    savefig("Figures/SVG/Comb - Count of Objects By Debris Status.svg")
end

function calculateBandedCollisionFrequency()
    # Calculate spatial density of all bands with all objects
    # Calculate collision frequency by band

    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    ALL_altitudes = getTleAltitude.(tles)
    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    filter!((x) -> x .< 2500 && x .> 100, ALL_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, DEB_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, NOD_altitudes)


    sliceWidth = 3
    bins = (0:sliceWidth:2500)

    ALL_binned = fit(Histogram, ALL_altitudes, bins, closed = :right)
    DEB_binned = fit(Histogram, DEB_altitudes, bins, closed = :right)
    NOD_binned = fit(Histogram, NOD_altitudes, bins, closed = :right)

    bins_arr = collect(bins)
    popat!(bins_arr, 1)

    ALL = ALL_binned.weights

    DEB_ADJ = DEB_binned.weights
    NONF = NOD_binned.weights .* 0.35 #Non functional satellites
    NOD_ADJ = NOD_binned.weights .* (1 - 0.35) #Functional satellites


    ALL_DENSITY = getSphereDensity.(bins_arr, ALL)
    DEB_ADJ_DENSITY = getSphereDensity.(bins_arr, DEB_ADJ)
    NONF_ADJ_DENSITY = getSphereDensity.(bins_arr, NONF)
    NOD_ADJ_DENSITY = getSphereDensity.(bins_arr, NOD_ADJ)

    Acc_DEB = 0.5 # Estimate
    Acc_SAT = 4 #Kessler
    Vs = 7 #Kessler

    # CF_ALL = ALL_DENSITY .^ 2 .* Acc_SAT .* Vs .* 

    # @show sum(CF_ALL)

    CF_ALL = 0.5 .* ALL_DENSITY .^ 2 .* Acc_SAT .* Vs .* ((4 / 3 * π) .* ((bins_arr .^ 3) .- (bins_arr .- sliceWidth) .^ 3))

    @show sum(CF_ALL)

    # @show sum(CF_ALL) / (4 / 3 * π * (2500^3))

    plot(bins_arr, CF_ALL, yscale = :log10, ylims = (10E-9, 10E-1), color = "darkred")

    title!("Rate of Collisions Per Year by Altitude")
    xlabel!(L"Altitude ($km$)")
    # xlims!(0, 1)
    ylabel!("Rate of Collisions")

    png("Figures/PNG/Comb - Rate of Collisions by Altitude")
    savefig("Figures/SVG/Comb - Rate of Collisions by Altitude.svg")
end

calculateBandedCollisionFrequency()
# Probability of Debris - Nonf collision, by altitude band
function calculateBandedCollisionFrequencyDebToNonf()
    # Calculate spatial density of all bands with all objects
    # Calculate collision frequency by band

    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    ALL_altitudes = getTleAltitude.(tles)
    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    filter!((x) -> x .< 2500 && x .> 100, ALL_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, DEB_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, NOD_altitudes)


    sliceWidth = 3
    bins = (0:sliceWidth:2500)

    ALL_binned = fit(Histogram, ALL_altitudes, bins, closed = :right)
    DEB_binned = fit(Histogram, DEB_altitudes, bins, closed = :right)
    NOD_binned = fit(Histogram, NOD_altitudes, bins, closed = :right)

    bins_arr = collect(bins)
    popat!(bins_arr, 1)

    ALL = ALL_binned.weights

    DEB_ADJ = DEB_binned.weights
    NONF = NOD_binned.weights .* 0.35 #Non functional satellites
    NOD_ADJ = NOD_binned.weights .* (1 - 0.35) #Functional satellites


    ALL_DENSITY = getSphereDensity.(bins_arr, ALL)
    DEB_ADJ_DENSITY = getSphereDensity.(bins_arr, DEB_ADJ)
    NONF_ADJ_DENSITY = getSphereDensity.(bins_arr, NONF)
    NOD_ADJ_DENSITY = getSphereDensity.(bins_arr, NOD_ADJ)

    Acc_DEB = 0.5 # Estimate
    Acc_SAT = 4 #Kessler
    Vs = 7 #Kessler

    # CF_ALL = ALL_DENSITY .^ 2 .* Acc_SAT .* Vs .* 

    # @show sum(CF_ALL)

    CF_DEB = 0.5 .* DEB_ADJ_DENSITY .^ 2 .* Acc_SAT .* Vs .* ((4 / 3 * π) .* ((bins_arr .^ 3) .- (bins_arr .- sliceWidth) .^ 3))

    CF_DEB ./ DEB_ADJ .* NONF # Scale to be collisions between nonf and deb rather than between deb

    @show sum(CF_DEB)

    # @show sum(CF_ALL) / (4 / 3 * π * (2500^3))

    plot(bins_arr, CF_DEB, yscale = :log10, ylims = (10E-9, 10E-1), legend = false, color = "darkred")

    title!("Rate of Collisions Per Year by Altitude, \n Debris to Nonfunctional Satellites")
    xlabel!(L"Altitude ($km$)")
    # xlims!(0, 1)
    ylabel!("Rate of Collisions")

    png("Figures/PNG/DEBTONONF - Rate of Collisions by Altitude")
    savefig("Figures/SVG/DEBTONONF - Rate of Collisions by Altitude.svg")
end


calculateBandedCollisionFrequencyDebToNonf()

#Probability of Nonf - Nonf collision, by altitude band
function calculateBandedCollisionFrequencyNonftoNonf()
    # Calculate spatial density of all bands with all objects
    # Calculate collision frequency by band

    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    ALL_altitudes = getTleAltitude.(tles)
    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    filter!((x) -> x .< 2500 && x .> 100, ALL_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, DEB_altitudes)
    filter!((x) -> x .< 2500 && x .> 100, NOD_altitudes)


    sliceWidth = 3
    bins = (0:sliceWidth:2500)

    ALL_binned = fit(Histogram, ALL_altitudes, bins, closed = :right)
    DEB_binned = fit(Histogram, DEB_altitudes, bins, closed = :right)
    NOD_binned = fit(Histogram, NOD_altitudes, bins, closed = :right)

    bins_arr = collect(bins)
    popat!(bins_arr, 1)

    ALL = ALL_binned.weights

    DEB_ADJ = DEB_binned.weights
    NONF = NOD_binned.weights .* 0.35 #Non functional satellites
    NOD_ADJ = NOD_binned.weights .* (1 - 0.35) #Functional satellites


    ALL_DENSITY = getSphereDensity.(bins_arr, ALL)
    DEB_ADJ_DENSITY = getSphereDensity.(bins_arr, DEB_ADJ)
    NONF_ADJ_DENSITY = getSphereDensity.(bins_arr, NONF)
    NOD_ADJ_DENSITY = getSphereDensity.(bins_arr, NOD_ADJ)

    Acc_DEB = 0.5 # Estimate
    Acc_SAT = 4 #Kessler
    Vs = 7 #Kessler

    CF_NONF = 0.5 .* NONF_ADJ_DENSITY .^ 2 .* Acc_SAT .* Vs .* ((4 / 3 * π) .* ((bins_arr .^ 3) .- (bins_arr .- sliceWidth) .^ 3))

    @show sum(CF_NONF)

    plot(bins_arr, CF_NONF, yscale = :log10, ylims = (10E-3, 10E2), legend = false, color = "darkred")

    title!("Rate of Collisions Per Year by Altitude, \n Non Functional Satellites")
    xlabel!(L"Altitude ($km$)")
    # xlims!(0, 1)
    ylabel!("Rate of Collisions")

    png("Figures/PNG/NONF - Rate of Collisions by Altitude")
    savefig("Figures/SVG/NONF - Rate of Collisions by Altitude.svg")
end

calculateBandedCollisionFrequencyNonftoNonf()

function doCalculateRunawayThreshold()
    tles = read_tle("TLEData.txt")

    ALL_altitudes = getTleAltitude.(tles)

    filter!((x) -> x .< 2500 && x .> 100, ALL_altitudes)

    bins = collect(0:50:2500)

    ALL_binned = fit(Histogram, ALL_altitudes, bins, closed = :right)

    CUMULATIVE = cumsum(reverse(ALL_binned.weights))

    popat!(bins, 1)


    Re = 6371 # radius of the earth, km
    G = 6.673 * 10^-11 # Gravitational constant
    Me = 5.972 * 10^24 # Mass of earth, Kg
    Cd = 2.2
    W = 1.5

    mA = 125 # Average mass over area - given by Kessler
    V = 7.5 # Avg. Relative veloctiy, km/sec, given by Kessler
    σf = 14

    # rNi_arr = (4 .* π .* (Re .+ bins) .^ 3 .* .√((G + Me) ./ (Re .+ bins)) .* expatmosphere.(bins .* 1000) .* Cd) ./ (CUMULATIVE .* W .* mA .* V .* σf)
    rNi_arr = (4 .* π .* (Re .+ bins) .^ 3 .* .√((G + Me) ./ (Re .+ bins)) .* 10^-11 .* Cd) ./ (CUMULATIVE .* W .* mA .* V .* σf)

    @show rNi_arr

    plot(bins, [rNi_arr ALL_binned.weights], labels = ["Critical Threshold" "Count of Objects"], color = ["darkred" "blue"])

    title!("Debris Population and Critical Threshold \n by Altitude")
    xlabel!(L"Altitude ($km$)")
    xlims!(100, 2500)
    ylabel!("Number of Debris Objects")

    png("Figures/PNG/COMB - Critical Threshold")
    savefig("Figures/SVG/COMB - Critical Threshold.svg")
end


doCalculateRunawayThreshold()

function doDebrisRankingModel()
    tles = read_tle("TLEData.txt")

    getTleAltitude.(tles)
end

function buildAllFigures()
    global eop = get_iers_eop()

    # The DCM (Direction Cosine Matrix) that rotates TEME into alignment with ITRF
    global D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), DatetoJD(2022, 1, 1, 0, 0, 0), eop)


    doInclinationAnalysis()

    doDebrisOnlyInclinationAnalysis()

    doDebrisAltitudeBinAnalysis()

    doDebrisSpacialDensityAnalysis()

    doProportionDebrisByAltitude()

    doCountDebrisNonDebrisComparisonByAltitude()

    calculateBandedCollisionFrequency()

    calculateBandedCollisionFrequencyDebToNonf()

    calculateBandedCollisionFrequencyNonftoNonf()
end

buildAllFigures()
