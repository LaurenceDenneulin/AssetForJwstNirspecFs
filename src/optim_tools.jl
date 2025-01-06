       
function cost_psf(x::AbstractArray{T,1},w::AbstractArray{T,1}, σ::T,  ρ::T) where {T<:AbstractFloat}
    @assert size(x) == size(w)
    # norm = maximum(x[max(round(Int64,ρ)-5,1):min(round(Int64,ρ)+5,length(x))])
    b = w .!= 0
    norm = maximum(b .* x)
    h=GaussianPSF(σ)
    f=zero(T)
    for k=1:length(x)
        r=x[k]/norm - h(k - ρ,0)
        wr=w[k]*r
        f += r*wr
    end
    return f
end

#=
function estimate_psf_parameters(d, w, pos)
   m,n = size(d)
    σinit = 1.25
    ρinit, kf, kl=slit2cam(d, pos)
    psf_param=zeros(2 ,n)
    par=[σinit, ρinit]
    for k=1:n
        if sum(d[:,k]) !=0
            par .=[σinit, min(par[2],ρinit)]
            par .= bobyqa(x->cost_psf(d[:,k], w[:,k], x[1], x[2]), par, rhobeg=1., rhoend=1e-8)[1]
            psf_param[:,k] .= par
        end
    end
return psf_param
end
=#
function estimate_psf_parameters(d::AbstractArray{T,2}, 
    w::AbstractArray{T,2}, 
    lambda::AbstractArray{T,2}, 
    pos::T) where {T<:AbstractFloat}

    psf_param=zeros(2, size(d,2))
    cent_init = findmax(mean(w .* d, dims=2))[2]
    prms_init = [.75, cent_init[1]]
    for i in axes(lambda,2)
        # select the lambda value at the initial center of the slit for each column
        l = cent_init[1]
        # iso-lambda indices
        iso_lam = findmin(x -> abs.(x .- lambda[l,i]), lambda; dims=2)[2]
        # data profile for the iso-lambda
        #FIXME: is there a better way than allocating a new vector for each
        #lambda? Maybe a view?
        r_prof = [d[n[1],n[2]] for n in iso_lam][:]
        if sum(r_prof) != 0. # the profile contains a PSF slice
            w_prof = [w[n[1],n[2]] for n in iso_lam][:]
            prms = bobyqa(x -> cost_psf(r_prof, w_prof, x[1], x[2]), 
                          prms_init, rhobeg=1.5, rhoend=1e-4)[1]
            psf_param[:,i] = prms
        end
    end
    return psf_param
end
#

function fit_polynomial!(d, w, fit, x; order=2)
    xx=[w .* x.^k for k=order:-1:0]
    A=zeros(order+1, order+1)
    b=zeros(order+1)

    for k=1:(order+1)
       for l=1:(order+1)
        A[k,l]=xx[k]'*xx[l]
       end
       b[k] = xx[k]'*d
    end         
    coef = A\b

    xx=[x.^k for k=order:-1:0]
    for k=1:(order+1)
        fit.+= coef[k]*xx[k]     
    end
    return coef
end

function fit_polynomial!(d, w, fit; order=2)
   n=length(d)
   return fit_polynomial!(d, w, fit, collect(1:n); order=order)
end


function robust_fit_polynomial!(d, w, fit, x; order=2, tol=1e-6, threshold=1.5)
    fit_temp=copy(fit)
    test_tol=1.
    coef=zeros(order+1);
    while test_tol>tol
        fill!(fit,0.)
        coef.= fit_polynomial!(d, w, fit, x; order=order)   
        Χ²=(d - fit).^2
        w = w .*( Χ².<threshold)
        test_tol=sum(abs.(fit - fit_temp))
        fit_temp .=copy(fit)
    end        
    return coef
end

function robust_fit_polynomial!(d, w, fit; order=2, threshold=1.5)
   n=length(d)
   return robust_fit_polynomial!(d, w, fit, collect(1:n); order=order, threshold=threshold)
end

