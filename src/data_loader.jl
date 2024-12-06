"""
    data, wgt, bpm, lambda, pos =  load_data(filename)
    
    where the input 'filename' is a String, and the output:
    
    -'data' is the data map;
    
    -'wgt' is the weight map;
    
    -'bpm' is the defective pixel map with (1=not defective, 0=defective);
    
    -'lambda' is the wavelength map;
    
    -'pos' is the position of the object in the slit (-1.5 being the top et 1.5 the bottom).

"""
function read_data(filename::AbstractString)
    
    hdr = read(FitsHeader,filename,ext=2)
    pos_vect=hdr["SRCYPOS"].float

    data=Float64.(readfits(filename,ext=2))'
    data[isnan.(data)] .= 0.
    dq=Float64.(readfits(filename,ext=3))'
    dq[isnan.(dq)] .= 0.
    err=Float64.(readfits(filename,ext=4))'
    err[isnan.(err)] .= 0.
    lambda=Float64.(readfits(filename,ext=5))'
    lambda[isnan.(lambda)] .= 0.

    bpm = dq .<=6;
    wgt = bpm ./(err.^2) 
    wgt[isnan.(wgt)] .= 0.
    wgt[wgt .== Inf] .= 0.

return data, wgt, bpm, lambda, pos_vect
end


"""

    data, wgt, lambda, ρ, cent_init =  load_data(dir;...)
    
    where the inputs:

    -'dir' is the path to the directory where the '_cal' data are stored (one directory per filter).
    
    the optional inputs:
    -'order' is the order of the polynomial fitting for the geometrical_calibration;
    -'sky_sub=true' performs a background substraction,
    
    and the outputs:
    -'data' is a data cube;
    -'wgt' is the weights cube;
    -'lambda' is the wavlengths cube;
    -'ρ' is the spatial position cube (in pixel);
    -'cent_init' is the initial position of the object in the slit (in pixel).
    
"""
function load_data(dir::AbstractString;
                   order::Integer=2,
                   sky_sub=false,
                   save=false)
                   
    if save
       mkpath("save/")
    end               
    filenames=readdir(dir)
    filenames=filenames[contains.(filenames,"_cal")]     
    nfiles=length(filenames)          
    d = read_data(dir*filenames[1])[1];
    m,n = size(d)
    data=zeros(m,n,nfiles)
    wgt=zeros(m,n,nfiles)
    bpm_map=zeros(m,n,nfiles)
    lambda=zeros(m,n,nfiles)
    ρ=zeros(m,n,nfiles)
    pos=zeros(nfiles)
    for i=1:nfiles
         d, w, bpm, λ, p = read_data(dir*filenames[i])
         data[:,:,i] .=d
         wgt[:,:,i] .= w 
         bpm_map[:,:,i] .=bpm 
         lambda[:,:,i] .= λ 
         pos[i] = p       
         ρ[:,:,i] .= geometric_calibration(data[:,:,i], w.*bpm_map[:,:,i], pos[i]; order=order, save=save)[1]
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
    
    if save
       
        numpos=findfirst("_nrs",filenames[1])[1]-1
        extpos=findfirst(".",filenames[1])[1]
        filename_beg=filenames[1][1:numpos-1]
        filename_end=filenames[1][numpos+1:extpos-1]
        
        writefits!("save/"*filename_beg*"d"*filename_end*".fits",["TYPE" => "data"], data);
        writefits!("save/"*filename_beg*"w"*filename_end*".fits",["TYPE" => "weights"],wgt);
        writefits!("save/"*filename_beg*"l"*filename_end*".fits",["TYPE" => "wavelength"], lambda);
        writefits!("save/"*filename_beg*"x"*filename_end*".fits",["TYPE" => "spatial position"],ρ);
    end

    return data, wgt, lambda, ρ, cent_init 
end

