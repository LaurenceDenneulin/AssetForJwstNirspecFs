

"""
    data, wgt, bpm, lambda, pos =  load_data(filename)
    
    where the input 'filename' is a String, and the output:
    
    -'data' is the data map;
    
    -'wgt' is the weight map;
    
    -'bpm' is the defective pixel map with (1=not defective, 0=defective);
    
    -'lambda' is the wavelength map;
    
    -'pos' is the position of the object in the slit (-1.5 being the top et 1.5 the bottom).

"""
function load_data(filename::AbstractString)
    clear_vect=[]
    pos_vect=[]
    FitsIO(filename, "r") do io
       d = read(FitsImage,io;ext=2);
       push!(pos_vect, d["SRCYPOS"]);
       data = Float64.(convert(Array,d))
       data[isnan.(data)] .= 0. 
       push!(clear_vect, data)
        for k=3:5
           data = Float64.(read(FitsArray,io;ext=k));
           data[isnan.(data)] .= 0. 
           push!(clear_vect, data)
        end
    end

    data=Float64.(clear_vect[1]'[:,:]);
    m,n = size(data);

    bpm=Float64.(clear_vect[2]' .<=6);
    bpm .*=  Float64.(data .> 0.)
    err=Float64.(clear_vect[3]'[:,:]);
    lambda=Float64.(clear_vect[4]'[:,:]);
    wgt = bpm ./(err.^2) ;
    wgt[isnan.(wgt)] .= 0.;
    wgt[wgt .== Inf] .= 0.;

return data, wgt, bpm, lambda, pos_vect[1]
end


"""

    data, wgt, lambda, rho, cent_init =  load_data(filename_beg, filename_end, nfiles;...)
    
    where the inputs:

    -'filename_beg' is the first part of the String before the integration number;
    -'filename_end' is the last part of the String after the integration number;
    -'nfiles' is the amount of integrations,
    
    the optional inputs:
    -'order' is the order of the polynomial fitting for the geometrical_calibration;
    -'sky_sub=true' performs a background substraction,
    
    and the outputs:
    -'data' is a data cube;
    -'wgt' is the weights cube;
    -'lambda' is the wavlengths cube;
    -'rho' is the spatial position cube (in pixel);
    -'cent_init' is the initial position of the object in the slit (in pixel).
    
"""
function load_data(filename_beg::AbstractString, 
                   filename_end::AbstractString, 
                   nfiles::Integer; 
                   order::Integer=2,
                   sky_sub=false)
    d = load_data(filename_beg*"1"*filename_end)[1];
    m,n = size(d)
    data=zeros(m,n,nfiles)
    wgt=zeros(m,n,nfiles)
    bpm_map=zeros(m,n,nfiles)
    lambda=zeros(m,n,nfiles)
    rho=zeros(m,n,nfiles)
    pos=zeros(nfiles)
    for i=1:nfiles
         d, w, bpm, λ, p = load_data(filename_beg*"$i"*filename_end)
         data[:,:,i] .=d
         wgt[:,:,i] .= w 
         bpm_map[:,:,i] .=bpm 
         lambda[:,:,i] .= λ 
         pos[i] = p       
         rho[:,:,i] .= geometric_calibration(data[:,:,i], bpm_map[:,:,i], pos[i]; order=order)[1]
    end
    
    if sky_sub
        data_nobg = copy(data)
        wgt_nobg=copy(wgt)
    
        for i=1:nfiles
            ind=cat(collect(1:i-1),collect(i+1:nfiles),dims=1)
            data_nobg[:,:,i] .-= sum(data[:,:,ind],dims=3)/(nfiles-1)
            wgt_nobg[:,:,i] .= 1 ./ (1 ./wgt[:,:,i] + sum(1 ./wgt[:,:,ind], dims=3)/(nfiles-1)^2 )
        end
        data .= data_nobg
        wgt .= wgt_nobg
    end
    wgt[isnan.(wgt)] .= 0.;
    wgt[wgt .== Inf] .= 0.;
    cent_init=slit2cam(data, pos)
    return data, wgt, lambda, rho, cent_init 
end

