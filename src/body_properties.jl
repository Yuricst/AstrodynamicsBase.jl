"""
Wrapper for accessing gm information from de431
"""

"""
    get_gm(naifID)

Function returns GM value of body specified by NAIF ID.
For body names, refer to: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

# Arguments
    - `naifIDs::Int`: tuple containing strings of naif ID to use to extract GM values. Multiple naifIDs may be passed in a single function call.

# Returns
    - `Float64`: GM value

# Examples:
```julia-repl
julia> get_gm(399)
398600.435436096
```
"""
function get_gm(naifID::Int)
    # call gm values
    de431gm = get_gm_de431()
    bdyname = "BODY" * string(naifID) * "_GM"
    return de431gm[bdyname]
end


"""
    get_gm_de431()

Function returns tuple of gm values from de431 spice kernel

Args:
    None
Returns:
    (dict): dictionary with fields defined by "BODY" + <NAIF body ID> + "_GM"", which contains a tuple of the GM value of the corresponding body
"""
function get_gm_de431()
    de431 = Dict(
        "BODY1_GM" => 2.2031780000000021E+04,
        "BODY2_GM" => 3.2485859200000006E+05,
        "BODY3_GM" => 4.0350323550225981E+05,
        "BODY4_GM" => 4.2828375214000022E+04,
        "BODY5_GM" => 1.2671276480000021E+08,
        "BODY6_GM" => 3.7940585200000003E+07,
        "BODY7_GM" => 5.7945486000000080E+06,
        "BODY8_GM" => 6.8365271005800236E+06,
        "BODY9_GM" => 9.7700000000000068E+02,
        "BODY10_GM" => 1.3271244004193938E+11,
        "BODY199_GM" => 2.2031780000000021E+04,
        "BODY299_GM" => 3.2485859200000006E+05,
        "BODY399_GM" => 3.9860043543609598E+05,
        "BODY499_GM" => 4.282837362069909E+04,
        "BODY599_GM" => 1.266865349218008E+08,
        "BODY699_GM" => 3.793120749865224E+07,
        "BODY799_GM" => 5.793951322279009E+06,
        "BODY899_GM" => 6.835099502439672E+06,
        "BODY999_GM" => 8.696138177608748E+02,
        "BODY301_GM" => 4.9028000661637961E+03,
        "BODY401_GM" => 7.087546066894452E-04,
        "BODY402_GM" => 9.615569648120313E-05,
        "BODY501_GM" => 5.959916033410404E+03,
        "BODY502_GM" => 3.202738774922892E+03,
        "BODY503_GM" => 9.887834453334144E+03,
        "BODY504_GM" => 7.179289361397270E+03,
        "BODY505_GM" => 1.378480571202615E-01,
        "BODY601_GM" => 2.503522884661795E+00,
        "BODY602_GM" => 7.211292085479989E+00,
        "BODY603_GM" => 4.121117207701302E+01,
        "BODY604_GM" => 7.311635322923193E+01,
        "BODY605_GM" => 1.539422045545342E+02,
        "BODY606_GM" => 8.978138845307376E+03,
        "BODY607_GM" => 3.718791714191668E-01,
        "BODY608_GM" => 1.205134781724041E+02,
        "BODY609_GM" => 5.531110414633374E-01,
        "BODY610_GM" => 1.266231296945636E-01,
        "BODY611_GM" => 3.513977490568457E-02,
        "BODY615_GM" => 3.759718886965353E-04,
        "BODY616_GM" => 1.066368426666134E-02,
        "BODY617_GM" => 9.103768311054300E-03,
        "BODY701_GM" => 8.346344431770477E+01,
        "BODY702_GM" => 8.509338094489388E+01,
        "BODY703_GM" => 2.269437003741248E+02,
        "BODY704_GM" => 2.053234302535623E+02,
        "BODY705_GM" => 4.319516899232100E+00,
        "BODY801_GM" => 1.427598140725034E+03,
        "BODY901_GM" => 1.058799888601881E+02,
        "BODY902_GM" => 3.048175648169760E-03,
        "BODY903_GM" => 3.211039206155255E-03,
        "BODY904_GM" => 1.110040850536676E-03,
        "BODY2000001_GM" => 6.3130000000000003E+01,
        "BODY2000002_GM" => 1.3730000000000000E+01,
        "BODY2000003_GM" => 1.8200000000000001E+00,
        "BODY2000004_GM" => 1.7289999999999999E+01,
        "BODY2000006_GM" => 9.3000000000000005E-01,
        "BODY2000007_GM" => 8.5999999999999999E-01,
        "BODY2000010_GM" => 5.7800000000000002E+00,
        "BODY2000015_GM" => 2.1000000000000001E+00,
        "BODY2000016_GM" => 1.8100000000000001E+00,
        "BODY2000029_GM" => 8.5999999999999999E-01,
        "BODY2000052_GM" => 1.5900000000000001E+00,
        "BODY2000065_GM" => 9.1000000000000003E-01,
        "BODY2000087_GM" => 9.8999999999999999E-01,
        "BODY2000088_GM" => 1.0200000000000000E+00,
        "BODY2000433_GM" => 4.463E-4,
        "BODY2000511_GM" => 2.2599999999999998E+00,
        "BODY2000704_GM" => 2.1899999999999999E+00,
    )
    return de431
end


"""
    get_body_radius(naifID::Int)

Function returns body radius from NAIF ID, in km
"""
function get_body_radius(naifID::Int)
    # planet radius
    radius_dict = Dict(
        "1" => 4879.0 / 2,
        "2" => 12104 / 2,
        "3" => 6378.0,
        "4" => 6792.0 / 2,
        "5" => 142984.0 / 2,
        "6" => 120536.0 / 2,
        "7" => 51118 / 2,
        "8" => 49528.0 / 2,
        "9" => 2370 / 2,
        "199" => 4879.0 / 2,
        "299" => 12104 / 2,
        "399" => 6378.0,
        "499" => 6792.0 / 2,
        "599" => 142984.0 / 2,
        "699" => 120536.0 / 2,
        "799" => 51118 / 2,
        "899" => 49528.0 / 2,
        "999" => 2370 / 2,

        # Jovian satellites
        "501" => 1821.5,   # Io
        "502" => 1560.8,   # Europa
        "503" => 2631.2,   # Ganymede
        "504" => 2410.3,   # Callisto
    )
    return radius_dict[string(naifID)]
end


"""
    get_body_period(naifID)

Function returns period from NAIF ID, in seconds.
Conversion: 1 year = 365.25 days = 31557600 s
Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/
"""
function get_body_period(naifID::Int)
    period_dict = Dict(
        "1" => 0.2408467,
        "2" => 0.61519726,
        "3" => 1.0000174,
        "4" => 1.8808158,
        "5" => 11.862615,
        "6" => 29.447498,
        "7" => 84.016846,
        "8" => 164.79132,
        "9" => 248.0208,
    )
    return period_dict[string(naifID)] * 31557600
end


"""
    get_body_period(naifID::Int, et_ref::Float64)

Function returns period from NAIF ID, in seconds.
Conversion: 1 year = 365.25 days = 31557600 s
Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/
"""
function get_body_period(naifID::Int, et_ref::Float64, center_body::Int)
    # get state-vector and compute SMA
    sv = spice_spkssb(et_ref, naifID, "ECLIPJ2000")
    if center_body == 10
        return get_period(sv, get_gm(center_body))
    else
        # shift origin
        sv_shifted = shift_eclipj2000_to_bodycenter(sv, center_body, et_ref)
        return get_period(sv_shifted, get_gm(center_body))
    end
end


"""
    get_body_sma(naifID)

Function returns SMA from NAIF ID, in km
Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/
"""
function get_body_sma(naifID::Int)
    au = 1.495978707e8
    sma_dict = Dict(
        "1" => 0.3870993,
        "2" => 0.723336,
        "3" => 1.000003,
        "4" => 1.52371,
        "5" => 5.2029,
        "6" => 9.537,
        "7" => 19.189,
        "8" => 30.0699,
        "9" => 39.4821,
    )
    return sma_dict[string(naifID)] * au
end


"""
    get_body_sma(naifID::Int, et_ref::Float64)

Function returns SMA from NAIF ID, in km
Data from: https://www.princeton.edu/~willman/planetary_systems/Sol/
"""
function get_body_sma(naifID::Int, et_ref::Float64, center_body::Int)
    # get state-vector and compute SMA
    sv = spice_spkssb(et_ref, naifID, "ECLIPJ2000")
    if center_body == 10
        return get_semiMajorAxis(sv, get_gm(center_body))
    else
        # shift origin
        sv_shifted = shift_eclipj2000_to_bodycenter(sv, center_body, et_ref)
        return get_semiMajorAxis(sv_shifted, get_gm(center_body))
    end
end


"""
    get_body_soi(naifID::Int, center_body::Int)

Compute Sphere of Influence (SOI) of body
"""
function get_body_soi(naifID::Int, center_body::Int)
    sma = get_body_sma(naifID)
    return sma * (get_gm(naifID) / get_gm(center_body))^(2 / 5)
end


"""
    get_canonical_param(
        center::Int = 10,
        ref_body::Int = 399,
        et_ref::Real = 6.311088691839073e8,
    )

Get canonical length, time, and velocity scales
"""
function get_canonical_param(
    center::Int = 10,
    ref_body::Int = 399,
    et_ref::Real = 6.311088691839073e8,
)
    if center == 10 && ref_body == 399
        mu_sun = 1.3271244004193938E+11
        au = 1.495978707e8
        lstar = au
        vstar = sqrt(mu_sun / au)
        tstar = au / vstar
    else
        mu_center = get_gm(center)
        sma = get_body_sma(ref_body, et_ref, center)
        lstar = sma
        vstar = sqrt(mu_center / sma)
        tstar = sma / vstar
    end
    return lstar, tstar, vstar
end

