using SatelliteToolbox
using Random
using Plots
using StatsBase
using LaTeXStrings
using DelimitedFiles

global eop = get_iers_eop()
# The DCM (Direction Cosine Matrix) that rotates TEME into alignment with ITRF
global D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), DatetoJD(2022, 1, 1, 0, 0, 0), eop)
gr()

function buildModel()
    tles = read_tle("TLEData.txt")

    tles_DEBRIS = filter(e -> occursin("DEB", e.name), tles) #Debris
    tles_NOD = filter(e -> !occursin("DEB", e.name), tles)  # Not Debris

    ALL_altitudes = getTleAltitude.(tles)
    DEB_altitudes = getTleAltitude.(tles_DEBRIS)
    NOD_altitudes = getTleAltitude.(tles_NOD)

    ALL_Comb = [tles ALL_altitudes]
    ALL_Comb = ALL_Comb[(ALL_Comb[:, 2].<2500).&(ALL_Comb[:, 2].>100), :]

    tles = ALL_Comb[:, 1]
    ALL_altitudes = ALL_Comb[:, 2]

    DEB_Comb = [tles_DEBRIS DEB_altitudes]
    DEB_Comb = DEB_Comb[(DEB_Comb[:, 2].<2500).&(DEB_Comb[:, 2].>100), :]

    tles_DEBRIS = DEB_Comb[:, 1]
    DEB_altitudes = DEB_Comb[:, 2]

    NOD_Comb = [tles_NOD NOD_altitudes]
    NOD_Comb = NOD_Comb[(NOD_Comb[:, 2].<2500).&(NOD_Comb[:, 2].>100), :]

    tles_NOD = NOD_Comb[:, 1]
    NOD_altitudes = NOD_Comb[:, 2]

    TLEResults = [analyzeTLE.(tles_NOD, Ref(tles_DEBRIS), Ref(tles_NOD), Ref(ALL_altitudes), Ref(DEB_altitudes), Ref(NOD_altitudes))][1]

    TLEResults = [getproperty.(tles_NOD, :name) mapreduce(permutedims, vcat, results)]

    @show TLEResults

    TLERankings = [TLEResults[:, 1] ordinalrank(TLEResults[:, 2]) ordinalrank(TLEResults[:, 3])]
    TLERankings = [TLERankings (TLERankings[:, 2] .+ TLERankings[:, 3])]

    TLERankings = sortslices(TLERankings, by = (x) -> x[4], dims = 1)

    TLERankings[:, 1] .= chop.(TLERankings[:, 1], head = 2, tail = 0)

    TLERankings = [["Satellite" "P Rank" "Alpha Rank" "Total Score"]; TLERankings]

    @show TLERankings
    writedlm("ModelOutput.csv", TLERankings, ",")
    writedlm("Top100ModelOutput.csv", TLERankings[1:102, :], ",")
end

# Helper Functions

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

function getSpatialDensity(radius, count, ΔA)
    outerVolume = 4 / 3 * π * radius^3
    innerVolume = 4 / 3 * π * (radius - ΔA)^3

    sliceVolume = outerVolume - innerVolume

    return count / sliceVolume
end

function getRunawayThreshold(altitude, ALL_altitudes)
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

    return rNi_arr[Int(fld(altitude, 50))]

end

#tles_NOD, Ref(tles_DEBRIS), Ref(tles_NOD), Ref(ALL_altitudes), Ref(DEB_altitudes), Ref(NOD_altitudes)
function analyzeTLE(tle, tles_DEBRIS, tles_NOD, ALL_altitudes, DEB_altitudes, NOD_altitudes)

    tle_alt = getTleAltitude(tle)
    alt_band = 5 # 10 km band, 5 on either side

    Vs = 7 #km / s
    Ac = 4 #m^2

    Qd = length(filter((x) -> (x .> (tle_alt - alt_band) && x .< (tle_alt + alt_band)), DEB_altitudes))
    Qs = length(filter((x) -> (x .> (tle_alt - alt_band) && x .< (tle_alt + alt_band)), NOD_altitudes))

    Sd = getSpatialDensity(tle_alt + alt_band, Qd, alt_band * 2)
    Ss = getSpatialDensity(tle_alt + alt_band, Qs, alt_band * 2) * 0.35
    Ps = (Sd + Ss) * Vs * Ac

    ttl_objs = Qd + Qs
    σd = Qd / ttl_objs
    σs = Qs / ttl_objs

    Xd = 400
    Xs = 1000

    Di = Xd * σd + Xs * σs

    Dp = Di * Ps

    ΔS = getSpatialDensity(tle_alt + alt_band, Qd + Dp, alt_band * 2) - Sd

    outerVolume = 4 / 3 * π * (tle_alt + alt_band)^3
    innerVolume = 4 / 3 * π * (tle_alt + alt_band - 10)^3

    ΔV = outerVolume - innerVolume

    ΔPf = 1 / 2 * ΔS^2 * Vs * Ac * ΔV

    rNi = getRunawayThreshold(tle_alt, ALL_altitudes)

    α = Dp / rNi

    return [ΔPf, α]
end



buildModel()


