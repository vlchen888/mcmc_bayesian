function data = moondata(d,N,sigma)
    %theta = unifrnd(0,2*pi,N,1);
    theta = (0:2*pi/N:2*pi-1/N)';
    data_clean = [cos(theta),sin(theta)] + (theta>pi)*[2,0];    
    data = [data_clean,zeros(N,d-2)] + sigma*normrnd(0,1,N,d);
end