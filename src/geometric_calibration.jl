

"""

"""
function geometric_calibration(d, bpm, pos; order=2)
    m,n = size(d)    
    psf_param=estimate_psf_parameters(d, bpm, pos)
    if order == 0
        x=[median(psf_param[2,:])]
    else        
        valid= Float64.(psf_param[2,:] .!=0)
        rho_shift = zeros(n)
        fit_polynomial!(psf_param[2,:], valid, rho_shift; order=order)
    end
    rho_shift .-=maximum(rho_shift);
    
    rho= Float64.(repeat(1:m, outer=(1,n)));      
    for k=1:n
        rho[:,k] .-= rho_shift[k]
    end 
    writedlm("rho_$pos.txt", rho_shift)
    return rho, rho_shift
end


"""

"""
function extract_λ_val(ρ_map, λ_map, ρ0)
    @assert size(ρ_map)==size(λ_map)
    m,n=size(ρ_map)
    λ = zeros(n)
    for k=1:n
        ρrange=range(minimum(ρ_map[:,k]); length=m, stop=maximum(ρ_map[:,k]))
        Iλ = SparseInterpolator(ker, [ρ0], ρrange)
        λ[k] = (Iλ*λ_map[:,k])[1]   
    end
    return λ
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
    return ((-pos .+0.5)*(ue-le-1) .+le .+1), kf,kl
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

