       
function cost_psf(x::AbstractArray{T,1},w::AbstractArray{T,1}, σ::T,  ρ::T) where {T<:AbstractFloat}
    @assert size(x) == size(w)
    norm = maximum(x[(max(round(Int64,ρ)-5,1):min(max(1,round(Int64,ρ))+5,length(x)))])
    h=GaussianPSF(σ)
    f=zero(T)
    for k=1:length(x)
        r=x[k]/norm - h(k - ρ,0)
        wr=w[k]*r
        f += r*wr
    end
    return f
end


function estimate_psf_parameters(d, w, pos)
   m,n = size(d)
    σinit = 0.5
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

function estimate_psf_parameters(d::AbstractArray{T,2}, 
    w::AbstractArray{T,2}, 
    lambda::AbstractArray{T,2}, 
    pos::T) where {T<:AbstractFloat}

    bpm = (w .!= 0)
    psf_param=zeros(2, size(d,2))
    cent_init = slit2cam(d, pos)
    lambda_cent = findmax(mean(bpm .* d, dims=2))[2]
    prms_init = [0.75, cent_init[1]]
    xl=[0., 1.]
    xu=[3., size(d,1)]
    for i in axes(lambda,2)
        # select the lambda value at the initial center of the slit for each column
        l = lambda_cent[1]
        # iso-lambda indices
        iso_lam = findmin(abs.(lambda .- lambda[l,i]); dims=2)[2]
        # data profile for the iso-lambda
        r_prof = bpm[iso_lam] .* d[iso_lam]
        if sum(r_prof) != 0. # the profile contains a PSF slice
            w_prof = w[iso_lam]#[w[n[1],n[2]] for n in iso_lam][:]
            max_pos = findmax(r_prof.*w_prof)[2][1];
            prms_init .= [0.75, max_pos]
            #=
            if abs(prms_init[2] - max_pos) < 1.5
                prms_init .= [0.75, max_pos]
            elseif (i>2) && (abs(prms_init[2] - psf_param[2,i-2]) < 0.5)
                prms_init .= [0.75, prms_init[2]]
            else
                prms_init .= [0.75, cent_init[1]]
            end
            =#
            prms = bobyqa(x -> cost_psf(r_prof[:], w_prof[:], x[1], x[2]),
                          prms_init, xl=xl,xu=xu, rhobeg=.5, rhoend=1e-8)[1]
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

