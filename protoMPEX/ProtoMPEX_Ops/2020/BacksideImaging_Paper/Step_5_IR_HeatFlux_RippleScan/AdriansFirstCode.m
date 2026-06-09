% Adrian first computer code!!!!

clear all
clc

prompt = 'Hello Adrian, will you shower tonight? [Y/N] ';
str = input(prompt,'s');

% if isempty(str)
%     str = 'Y';
% end

switch str
    case 'Y'
        disp('Yaay!!! no more stinky tooshie!!!')
        disp('I have a present for you!')
        
        prompt = 'Do you want to see the present? [Y/N] ';
        str2 = input(prompt,'s');
        
        switch str2
            case 'Y'
                clc
                figure
                penny
                view([0,90])
                title('Money for Adrian!!')
            case 'N'
                disp('Ok no worries! have a good night!!')           
        end
        
    case 'N'
        disp('mmmmm stinky tooshie!!')
        for ss = 1:12
            beep 
            pause(0.5)
        end
        
    otherwise
        disp('Say what!!!???')

end


