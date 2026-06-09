% #####################################################################
    % EXCEPTIONS =============================================================
    % 2017_06_16: The mach probe 10.5 V was swept with the isolation transformer in
    % order to get a the difference in the ion saturation current. this was not
    % recorded on the logbook. shots: 15168, 15170, 15171
    if sum(Shot(s) == [15168, 15170, 15171])
        rngMean{s} = [];
    end
    % 2017_07_18: In the first 6 shots of the day we had 2 shots that were not
    % in the high density mode but they were recorded as HMJ = 1. 
    if sum(Shot(s) == [15714, 15715])
        HMJ(s) = 0;
    end
    % Other shots that were not in the high density mode
    if sum(Shot(s) == [15083,15084,14656,13918])
        HMJ(s) = 0;
    end
    % This shot has unsually high mean Te due to a glitch that could not be
    % detected by the code
    if sum(Shot(s) == [14272,14365,13857,13984,14130,14009,14189])
        rngMean{s} = rngMean{s}(1:end-2);
    elseif sum(Shot(s) == [13910])
        rngMean{s} = rngMean{s}(2:end);
    elseif sum(Shot(s) == [12645])
        rngMean{s} = rngMean{s}(1:end-3);    
    end
    % 2017_07_18: we connected the IF probe with an arbitrary off-axis target
    % DLP and swept it as a DLP. the IV traces are not very good
    if sum(Shot(s) == [15723:15725])
        HMJ(s) = 0;
    end
    % 2017_01_12: DLP 9.5 on this day have significanlty different measurements
    % when compared to DLP 10.5. Te measured to be 6 eV while DLP 10.5 was
    % measured to be 1.5 eV
    if sum(Shot(s) == [12791,12793])
        HMJ(s) = 0;
    end
        % 2017_01_05: DLP 4.5. we used the non-isolated sweep box for this
        % measurement and it distorted the IV trace and giving unsually
        % high Te
    if sum(Shot(s) == [12414,12416,12415])
        HMJ(s) = 0;
    end
    % =========================================================================
    % #####################################################################