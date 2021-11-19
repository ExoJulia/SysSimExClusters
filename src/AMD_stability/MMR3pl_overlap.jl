#=
Date: 8 October 2020

Author: A. C. Petit

Functions giving spacing limits for multiplanetary systems based on the overlap of zeroth order 3 planets systems.
Based on A. C. Petit, G. Pichierri, M.B. Davies, A. Johansen 2020 (https://ui.adsabs.harvard.edu/abs/2020A%26A...641A.176P/abstract)

=#

"""
    MMRoverlapcriterion(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real,fudge::Real=1)

Tests if the three planet MMR overlap and returns a symbol (`:Stable` or `:MMR3ploverlap`) as well as the ratio of `δ/δov` of the spacing over the overlap spacing (larger than 1 means stable).

νij is the period ratio between planet i and j, ϵi is the planet i-to-star mass ratio, `fudge` corresponds to a factor that can be added to increase the amount of resonances, simulating the effect of additionnal planets in the system.
"""
function MMRoverlapcriterion(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real,fudge::Real=1)
    α12 = ν12^(2/3)
    α23 = ν23^(2/3)
    δ = genspacing(α12,α23)
    
    δov = fudge^.25*spacingoverlap3planets(ν12,ν23,ϵ1,ϵ2,ϵ3)
    if δ>δov
        return :Stable, δ/δov
    else
        return :MMR3ploverlap, δ/δov
end


"""
    survivaltime(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real,fudge::Real=1)

Compute the survival time of a system according to the 3 pl MMR overlap (Eq. 81 Petit et al. 2020)
νij is the period ratio between planet i and j, ϵi is the planet i-to-star mass ratio, `fudge` corresponds to a factor that can be added to increase the amount of resonances, simulating the effect of additionnal planets in the system. A good estimate is fudge=1 for 3, fudge=2 for five planets, or more.

Returns the average time, in units of the innermost period. If the resonance network is not overlapped, returns Inf
"""
function survivaltime(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real,fudge::Real=1)

    α12 = ν12^(2/3)
    α23 = ν23^(2/3)
    δ = genspacing(α12,α23)
    
    δov = fudge^.25*spacingoverlap3planets(ν12,ν23,ϵ1,ϵ2,ϵ3)

    if δ>=δov
        return +Inf
    else
        ratioδ = δ/δov
        η = resonancelocator(ν12,ν23) 
        
        ϵM = massfactor(η,α12,α23,ϵ1,ϵ2,ϵ3)
        facA = √(38/π) #

        prefactor = ϵM*ν12*facA*√(η*(1-η))/fudge^2 #Of the diffusion coefficient

        Tnorm = 2^1.5/9*(ratioδ)^6/(1-ratioδ^4)*10^√(-log(1-ratioδ^4))

        Tsurv = (3/2)^2/prefactor*Tnorm*3/32 #Deta=3/2 in units of δ
        return Tsurv
    end
end

"""
    spacingoverlap3planets(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real)

Determine the generalized spacing (defined by genspacing) were a triplet of planets can become unstable due to the overlap of three planet MMR. νij is the period ratio between planet i and j, ϵi is the planet i-to-star mass ratio.

Eq. (59) Petit et al. 2020

Returns the generlaized overlap spacing.
"""
function spacingoverlap3planets(ν12::Real,ν23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real)
    @assert (0<=ν12<=1)&&(0<=ν23<=1) "Period ratios ν12 = $ν12 and ν23 = $ν23 should be smaller than 1."

    η = resonancelocator(ν12,ν23) 
    α12 = ν12^(2/3)
    α23 = ν23^(2/3)
    ϵM = massfactor(η,α12,α23,ϵ1,ϵ2,ϵ3)
    Aov = 4*√(2*38/π)/3 #
    return (ϵM*Aov*(η*(1-η))^(3/2))^(1/4)
end

"""
    resonancelocator(ν12::Real,ν23::Real)

Resonance locator η (Eq. 21 Petit et al. 2020)
"""
resonancelocator(ν12::Real,ν23::Real) = ν12*(1-ν23)/(1-ν12*ν23)

"""
    massfactor(η::Real,α12::Real,α23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real)
Mass factor for 3 planets overlap (Eq. 53 Petit et al. 2020). η is the resonance locator, αij is the semi-major axis ratio, ϵi is the planet i-to-star mass ratio
"""
massfactor(η::Real,α12::Real,α23::Real,ϵ1::Real,ϵ2::Real,ϵ3::Real) = √(ϵ1*ϵ3+ϵ2*ϵ3*(η/α12)^2+ϵ2*ϵ1*((1-η)*α23)^2)

"""
    genspacing(α12::Real,α23::Real)

Generalized spacing (Eq. 45 Petit et al. 2020), αij is the semi-major axis ratio of planet i and j
"""
genspacing(α12::Real,α23::Real) = (1-α12)*(1-α23)/(2-α12-α23)

