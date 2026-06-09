%%%Interp function specifically for the Phase Detection
function [y0,x1,x2,y1,y2] = interp_pp(x,y,x0)
y0 = zeros(length(x0),1);
for ii = 1:length(x0)
    xnew = x0(ii);
    %%Make sure these are not end values
    if xnew~=x(1) && xnew~=x(end)
        %%Make sure we are less then the maximum value
        if xnew<max(x)
            x1 = min(x(x>=xnew));
        end
        if xnew>max(x)
            x1 = max(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Catch 0 and 360%%%%
        if length(y(x==x1))==1
            y1 = y(x==x1);
        end
        if length(y(x==x1))~=1
            y1 = y(1);
        end
        %%%%%%%%%%%%%%%%%%%%%
        %%Make sure we are greater then the manimum value
        if xnew>min(x)
            x2 = max(x(x<=xnew));
        end
        if xnew<min(x)
            x2 = min(x);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Catch 0 and 360%%%%
        if length(y(x==x2))==1
             y2 = y(x==x2);
        end
        if length(y(x==x2))~=1
            y2 = y(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        if x1~=x2
        ynew = y1+(xnew-x1).*(y2-y1)./(x2-x1);
        end
        if x1==x2
            ynew = y1;
        end
    end
    if xnew==x(1)
        ynew = y(1);
    end
    if xnew==x(end)
        ynew = y(end);
    end
    y0(ii) = ynew;
end
