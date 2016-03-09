
params = [1.0000;6.6423];
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1;
figure;
for aa=0.1:0.1:1.9
%    for bb=3:0.1:10
        plot(Phit,qform(Phit,[aa,bb]),':'); hold on;
%    end
end
plot(Phit,qform(Phit,params),'-k');
