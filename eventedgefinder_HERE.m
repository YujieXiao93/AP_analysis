function up_idx = eventedgefinder_HERE(data,threshold,up,Num)
if size(data,2)==1;data=data';end


switch up
    case 1
        A=find(data>threshold);
    case -1
        A=find(data<threshold);
        
end

if ~isempty(A)
    idx= find(diff(A)~=1);
    start_idx = [1 idx+1];
    end_idx= [idx numel(A)];
    up_idx= A([start_idx ; end_idx]');
else
    up_idx=[];
end

if Num>0&~isempty(up_idx)
    up_idx=up_idx(1:Num,:);
end

end

