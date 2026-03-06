function p=recompactUp(ph)

N=size(ph,1);
L=size(ph,2)
for n=1:N
    while ph(n,1)>3
        ph(n,:)=ph(n,:)-2*pi;
    end
    for i=2:L
        while ph(n,i)<ph(n,i-1)-2
            ph(n,i)=ph(n,i)+2*pi;
        end
    end
end

p=ph;