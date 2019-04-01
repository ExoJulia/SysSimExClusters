##### Function to use the photoevaporation boundary in Carrera et al 2018 to separate planets:

function photoevap_boundary_Carrera2018(R::Float64, P::Float64)
    #R is the planet radius in Earth radii, P is the period in days
    #This function returns 1 if the planet is above the boundary, and 0 if the planet is below the boundary as defined by Eq. 5 in Carrera et al 2018
    Rtrans = 2.6*P^(-0.1467)
    if R > Rtrans
        return 1
    elseif R < Rtrans
        return 0
    end
end

