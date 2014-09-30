function [ Mult ] = create_spectrogram_multiplier( M,N,type,rs )
%CREATE_SPECTROGRAM_MULTIPLIER Create a spectrogram multiplier. 
%   Usage: [ Mult ] = create_spectrogram_multiplier( M,N,type,real_signal )
%
%   Input parameters:
%         M     : Number of frequency channel.
%         N     : Number of shift in time: L/a.
%         type  : Type of multiplier (string).
%         rs    : Boolean for real signal.
%         M     : Number of frequency channel.
%   
%   Output parameters:
%         Mult  : Output mulitplier in a matrix form.
%
%   This function create a spectrogram multiplier compatible with the LTFAT
%   toolbox.
%
%   The different type of multiplier are:
%
%   * low_pass      : 1 for 1/3 of lowest frequency, 0 otherwise
%   * high_pass     : 1 for 1/3 of highest frequency, 0 otherwise
%   * band_pass     : 1 for 1/3 of middle frequency, 0 otherwise
%   * very_low_pass : 1 for 1/6 of lowest frequency, 0 otherwise
%   * LTFAT         : LTFAT logo 
%   * time_pass     : cut the signal in time
%   * full          : Full 1 matrix -- Identity operator
%   * rand          : Random matrix with function rand
%   
%   See the code for other multiplier


if rs
    M2=M/2+1;
    Mult=zeros(M2,N);

    if strcmp(type,'low_pass')
        Mult(1:round(M/6)+1,:)=1;
    elseif strcmp(type,'high_pass')
        Mult(round(M/3):M2,:)=1;
    elseif strcmp(type,'band_pass')
        Mult(round(M/3):M2,:)=1;
        Mult(1:round(M/6)+1,:)=1;
        Mult=1-Mult;
    elseif strcmp(type,'very_low_pass')
        Mult(1:round(M/24)+1,:)=1;
    elseif strcmp(type,'LTFAT')
        Mult=ltfattext; % size = 401 * 600
        Mult = imresize(Mult, [M2 N]);
    elseif strcmp(type,'time_pass')
        Mult(:,N-round(N/3):N)=1;
        Mult(:,1:round(N/3))=1;
        Mult=1-Mult;
    elseif strcmp(type,'time_sine')
        Mult=repmat(sin((1:N)/N*2*pi*2),M2,1);
    elseif strcmp(type,'freq_sine')
        Mult=repmat(sin([1:M2]'/M2*2*pi*2),1,N);
    elseif strcmp(type,'freq_gauss')
        Mult=repmat(exp(-abs(([1:M2]'-M2/2)/(M2))),1,N);
    elseif strcmp(type,'full')
        Mult=ones(size(Mult));
    elseif strcmp(type,'rand')
        Mult=rand(size(Mult));
    else
       error('Unknown type of multiplier! \n')
    end
else
    Mult=zeros(M,N);

    if strcmp(type,'low_pass')
        Mult([1:round(M/6),M-round(M/6):M],:)=1;
    elseif strcmp(type,'high_pass')
        Mult(round(M/3):M-round(M/3),:)=1;
    elseif strcmp(type,'band_pass')
        Mult([1:round(M/6),M-round(M/6):M],:)=1;
        Mult(round(M/3):M-round(M/3),:)=1;
        Mult=1-Mult;
    elseif strcmp(type,'very_low_pass')
        Mult([1:round(M/24),M-round(M/24):M],:)=1;
    elseif strcmp(type,'time_pass')
        Mult(:,N-round(N/3):N)=1;
        Mult(:,1:round(N/3))=1;
        Mult=1-Mult;
    elseif strcmp(type,'time_sine')
        Mult=repmat(sin((1:N)/N*2*pi*2),M,1);
    elseif strcmp(type,'freq_sine')
        Mult=repmat(sin([1:M]'/M*2*pi*2),1,N);
    elseif strcmp(type,'freq_gauss')
        Mult=repmat(exp(-abs([1:M]'-M/2)/(M)),1,N);
    elseif strcmp(type,'full')
        Mult=ones(size(Mult));       
     elseif strcmp(type,'rand')
        Mult=rand(size(Mult));     
   elseif strcmp(type,'LTFAT')
        Mult=ltfattext; % size = 401 * 600
        Mult = imresize([Mult; flipud(Mult)], [M N]);
    else
       error('Unknown type of multiplier! \n')
    end
    
end

end
