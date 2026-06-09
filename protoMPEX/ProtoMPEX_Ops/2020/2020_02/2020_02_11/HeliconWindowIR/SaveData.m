% fileName must a string array 
% variableNames must be a cell

disp('Saving data...')
tic
    try
        save(fileName,variableNames{:},'-append');
    catch
        save(fileName,variableNames{:})
    end
toc
disp('Save complete!!');
whos ('-file',fileName')