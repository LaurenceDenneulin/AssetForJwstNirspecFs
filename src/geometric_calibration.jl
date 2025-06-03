"""

"""

# function geometric_calibration(d, bpm, pos; order=3,threshold=1.5, save=false)
function geometric_calibration(d, w, lambda, pos; order=3,threshold=1.5, save=false)
    m,n = size(d)    
    # psf_param=estimate_psf_parameters(d, bpm, pos)
    psf_param=estimate_psf_parameters(d, w, lambda, pos)
    valid= (psf_param[2,:] .>0)
    valid[1:15] .= 0
    valid[end-15:end] .= 0
    psf_param[2,:] .= min.(m, psf_param[2,:])
    psf_param[2,:] .= max.(1., psf_param[2,:])
    psf_centers = zeros(n)
    robust_fit_polynomial!(psf_param[2,:], valid, psf_centers; order=order, threshold=threshold)
    #coef = fit_polynomial_law((1,order), (1.0:n)[valid], psf_param[2,valid])
    #rho_pol = PolynLaw{1,Float64,order}(coef)
    #psf_centers = rho_pol.(1.0:n)

    if save
        writedlm("save/fwhm_$pos.txt", psf_param[1,:])
        writedlm("save/rho_$pos.txt", [psf_param[2,:] psf_centers])
    end
    display(maximum(psf_centers))
    rho_shift = psf_centers .-maximum(psf_centers);
    #=
    rho= Float64.(repeat(1:m, outer=(1,n)));      
    for k=1:n
        rho[:,k] .-= rho_shift[k]
    end 
    =#
    rho = zeros(m,n)
    id_pix = zeros(m,n)
    rho_rng = Float64.(1:m)
    for i in axes(lambda,2)
        l = psf_centers[i]
        iso_lam = findmin(abs.(lambda .- lambda[round(Int64,l),i]); dims=2)[2]
        rho[iso_lam] = rho_rng .- rho_shift[i]
        id_pix[iso_lam] .= 1
    end
    mask = (lambda .!=0.)
    rho .*= mask
    # fill the empty pixels that are between iso-lams (this is NOT the way...)
    empty_pix = findall((id_pix .== 0.) .* mask)
    for p in empty_pix
        if p[1] > 1
            rho[p] = (rho[p[1]-1,p[2]] + rho[p[1]+1,p[2]])/2
        else
            rho[p] = rho[p[1]+1,p[2]]
        end
    end
    #
    return rho, rho_shift#, psf_centers
end


function geometric_calibration(d, bpm, pos; order=3,threshold=1.5, save=false)
    m,n = size(d)    
    psf_param=estimate_psf_parameters(d, bpm, pos)
      
    valid= Float64.(psf_param[2,:] .>0)
    psf_param[2,:] .= min.(m, psf_param[2,:])
    #psf_param[2,:] .= max.(0., psf_param[2,:])
    psf_centers = zeros(n)
    robust_fit_polynomial!(psf_param[2,:], valid, psf_centers; order=order, threshold=threshold)
    
    if save
        writedlm("save/fwhm_$pos.txt", psf_param[1,:])
        writedlm("save/rho_$pos.txt", [psf_param[2,:] psf_centers])
    end
    #rho_shift = psf_centers .-maximum(psf_centers);

    rho= Float64.(repeat(1:m, outer=(1,n)));      
    for k=1:n
        rho[:,k] .-= psf_centers[k]#rho_shift[k]
    end 
    return rho#, rho_shift#, psf_centers
end



"""

"""
function slit2cam(d::AbstractArray{T,2}, pos::T) where {T<:AbstractFloat}
    left_corner = findfirst(d.!=0)
    kr = left_corner[1]
    kf = left_corner[2] + 10
    le=findfirst(d[:,kf] .!=0.)
    kr = min(kr,le)
    kl=findlast(d[kr,:] .!=0.) - 10
    ue=findlast(d[:,kf] .!=0.)
    kf=findfirst(d[kr,:] .!=0.) + 10
    return ((-pos .+0.55)*(ue-le-1)./1.1 .+le), kf,kl
end

"""

"""
function slit2cam(d::AbstractArray{T,3}, pos::AbstractVector{T}) where {T<:AbstractFloat}
    @assert length(pos) ==size(d)[3]
    ρinit = similar(pos)
    n = size(d)[2]
    for k = 1:length(pos)
       ρ,l,r = slit2cam(d[:,:,k],pos[k])
       ρinit[k] = ρ
       n = min(n, r-l+20)
    end

    return  ρinit,n
end

