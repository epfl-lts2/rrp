function main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                                solving_method,param,paramplot)
%MAIN main function    
%   This function creates and solves the problem.
%
%   See the code for details.



%% Load the sound
s=load_sound(sound_name);
s=check_sound(s);
Ls=length(s);                       % Length of the signal
L=dgtlength(Ls,a,M);                % Length of the transform
N= L/a;                             % Number of time shifts
% wavwrite(s,fs,strcat('comparaison/original_',sound_name,'.wav'));

%% Creation of the frame
fprintf('-- Create the windows and the operators... ');
g = create_window(a,M,L ,window_type );
    
[~,GB]= gabframebounds(g,a,M);
if real_signal
    G = @(x) dgtreal(x,g,a,M);
    Gt= @(x) idgtreal(x,g,a,M)/GB;
else
    G = @(x) dgt(x,g,a,M);
    Gt= @(x) idgt(x,g,a,M)/GB;
end

fprintf(' Done\n');


%% Creation of the Spectrogram Multiplier
fprintf('-- Create the spectrogram muliplier... ');
A=create_spectrogram_multiplier( M,N,type_multiplier,real_signal );
fprintf(' Done\n');

%% Starting point
S=G(s).*A; % Simple application of the spectrogram multiplier

%% Apply the spectrogram muliplier and reconstruction
nsim=0;
% Launch all the simulation with all methods and alpha possible.
for ii=1:length(solving_method)
    if strcmp(solving_method{ii},'FGLA')
        for jj=1:length(alpha)
            fprintf(strcat('-- Reconstruction method: ',solving_method{ii}, ...
             '  with alpha = %g\n'),alpha(jj));

            nsim=nsim+1;

            [ap,Ap,info_reconstruct] = spectrogram_reconstruction( S,G,Gt, ...
                       solving_method{ii},alpha(jj),param); %#ok<ASGLU>
            fsnr=ssnr(S, Ap);
            % wavwrite(ap,fs,['comparaison/reconstruction_' sound_name ... 
            %                   '_' solving_method{ii} '.wav']);
            R(nsim,:)=info_reconstruct.ssnr; %#ok<AGROW>
            fprintf('   * The obtained ssnr is: %g\n', fsnr);
            fprintf('-- Reconstruction done \n');
        end
    else
        fprintf(strcat('-- Reconstruction method: ',solving_method{ii},'\n'));
        [ap,Ap,info_reconstruct] = spectrogram_reconstruction( S,G,Gt, ...
                            solving_method(ii),alpha(1),param); %#ok<ASGLU>
        fsnr=ssnr(S, Ap);
        nsim=nsim+1;

        %wavwrite(ap,fs,['comparaison/reconstruction_' sound_name '_' ...
        %                                   solving_method{ii} '.wav']);
        R(nsim,:)=info_reconstruct.ssnr; %#ok<AGROW>
        fprintf('   * The obtained ssnr is: %g\n', fsnr);
        fprintf('-- Reconstruction done \n');
    end

end


%% Display the result



fprintf('-- Display the results... ');

nsim=0;
% arrange the name
for ii=1:length(solving_method)
    if strcmp(solving_method(ii),'FGLA')
        for jj=1:length(alpha)
            nsim=nsim+1;
            text_fig(nsim)={strcat('FGLA: alpha=',num2str(alpha(jj)))}; %#ok<AGROW>
        end
    else
        nsim=nsim+1;
        text_fig(nsim)=solving_method(ii); %#ok<AGROW>
    end
end

cfig=figure;

set(cfig, 'Position', paramplot.position)
set(gcf,'PaperPositionMode','auto')
iterations=1:param(1).maxit;
semilogx(iterations,R);
title(sound_name,'FontSize',14);
legend(text_fig,'Location',paramplot.legendlocation);
xlabel('Iterations ','FontSize',12);
ylabel('ssnr (dB)','FontSize',12);
drawnow;
fprintf(' Done\n');
if paramplot.save
    filename=strcat(paramplot.pathfigure,num2str(length(alpha)),'_', ...
                type_multiplier,'_',sound_name,'_',window_type);
    print('-dpng','-zbuffer','-r300',[filename,'.png']);
    hgsave([filename,'.fig'])
end
end
