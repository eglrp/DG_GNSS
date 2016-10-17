function cdt = clock_bias(M,Ar,u,r)
    n = length(u);
    % Ar*x = u solution
    x_ur = M*u;
    hr = ones(n,1);
    x_ur_last = (1/n)*hr'*(u - Ar*x_ur);
    x_u = [x_ur; x_ur_last];
    
    % Ar*x = r solution
    x_rr = M*r;
    x_rr_last = (1/n)*hr'*(r - Ar*x_rr);
    x_r = [x_rr; x_rr_last];
    
    % clock bias calculation
    if x_r'*[r;0] < 0
        cdt = (x_u'*[r;0] + x_r'*[u;1] + sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));
    else
        cdt = (x_u'*[r;0] + x_r'*[u;1] - sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));
    end
end