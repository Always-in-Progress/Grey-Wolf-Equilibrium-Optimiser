function nx = decouple(ns)
    N = numel(ns) / 24;
    nx=zeros(24,N);
    for i=1:24
        nx(i,:)=ns(1,(N*(i-1)+1):(N*i))';
    end
end