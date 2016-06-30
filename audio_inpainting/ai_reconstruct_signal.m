function srec = ai_reconstruct_signal(G, shole, transitions, param)

time = G.time;
ratio = G.ratio;
fs = G.fs;

srec = cell(numel(transitions),1);
for ii = 1:numel(transitions)
    transition = transitions{ii};
    if param.keeplength

        sdiff =  time(transition(2,2))-time(transition(1,1))-...
            time(transition(2,1)) + time(transition(1,2));
        srec{ii} = shole;
        param.fs = fs;
        shift = ai_fineadjust(srec{ii},...
                           [time(transition(1,1)), time(transition(1,2));...
                           time(transition(2,1))+sdiff,time(transition(2,2))],...
                           ratio,param);
        switch param.melt
            case 0
                srec{ii}( time(transition(1,1)):time(transition(2,2)) ) = ...
                     shole( shift+(time(transition(1,2)):(time(transition(2,1))+sdiff)));
            case 1
                len = length(time(transition(1,1)):time(transition(2,2)));
                win = ones(len,1);
                trans = floor(fs/8);
                win([len-trans+1:len,1:trans]) = firwin('hann',2*trans);

                srec{ii}( time(transition(1,1)):time(transition(2,2)) ) = ...
                    win.*shole( shift+(time(transition(1,2)):(time(transition(2,1))+sdiff))) + ...
                    (1-win).*shole( time(transition(1,1)):time(transition(2,2)) );
            case 2
                srec{ii} = ai_tf_reconst(shole,...
                    [time(transition(1,1)), shift+time(transition(1,2));...
                    shift+time(transition(2,1))+sdiff,time(transition(2,2))],...
                    ratio,param);
            otherwise
                error('Unknown cross-fading method');
        end    
     
    else
        param.fs = fs;
        [shift,shiftr] = ai_fineadjust(shole,...
                           [time(transition(1,1)), time(transition(1,2)); ...
                           time(transition(2,1)),time(transition(2,2))],...
                           ratio,param);
        switch param.melt
            case 0
                srec{ii} = ...
                [shole(1:time(transition(1,1))); ...
                    shole(shift+(time(transition(1,2)):time(transition(2,1))));...
                    shole( time(transition(2,2)):end)];
            case 1
                trans = floor(fs/8);
                len = length(time(transition(1,2)):time(transition(2,1)))+2*trans;
                win = ones(len,1);
                trans = floor(fs/8);
                win([len-trans+1:len,1:trans]) = firwin('hann',2*trans);
                temp = [shole(time(transition(1,1))+(-trans:-1));zeros(len-2*trans,1);...
                    shole(time(transition(2,2))+(1:trans))].*(1-win);
                temp = temp+win.*shole(shift+(time(transition(1,2)) ...
                    -trans:time(transition(2,1))+trans));
                srec{ii} = ...
                [shole(1:time(transition(1,1))-trans); ...
                    temp;shole(trans+time(transition(2,2)):end)];
            case 2
                srec{ii} = ai_tf_reconst(shole,...
                    [time(transition(1,1)), shift+time(transition(1,2));...
                    time(transition(2,1)),time(transition(2,2))-shiftr],...
                    ratio,param);
            otherwise
                error('Unknown cross-fading method');
        end
    end
end

end