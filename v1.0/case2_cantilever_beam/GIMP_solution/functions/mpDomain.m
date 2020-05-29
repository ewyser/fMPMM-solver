function[mpD] = mpDomain(l,mpD);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ip=find(mpD.l(:,1)> l.*mpD.l0(:,1));
jp=find(mpD.l(:,1)<=l.*mpD.l0(:,1));
np=2*length(ip)+length(jp);

if(isempty(ip)<1)
    lpx = 0.5.*mpD.l(ip,1).*[-1 1];
    lpy = mpD.l(ip,2).*[1 1];
    
    mpD.x(1:np,1)  = [mpD.x(jp,1)   ;reshape(mpD.x(ip,1)+lpx,2*size(ip,1),1)     ];
    mpD.x(1:np,2)  = [mpD.x(jp,2)   ;reshape([mpD.x(ip,2) mpD.x(ip,2)],2*size(ip,1),1)     ];
    mpD.v(1:np,1) = [mpD.v(jp,1) ; reshape(repmat(mpD.v(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.v(1:np,2) = [mpD.v(jp,2) ; reshape(repmat(mpD.v(ip,2),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,1) = [mpD.p(jp,1) ; reshape(repmat(mpD.p(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,2) = [mpD.p(jp,2) ; reshape(repmat(mpD.p(ip,2),1,2) ,2*size(ip,1),1)];
    
    mpD.n0 = [mpD.n0(jp)  ; reshape(repmat(mpD.n0(ip),1,2),2*size(ip,1),1)];
    mpD.l  = [mpD.l(jp,:) ; [abs(lpx(:)) abs(lpy(:))]                  ];
    mpD.l0 = [mpD.l0(jp,:); [abs(lpx(:)) abs(lpy(:))]                  ];
    mpD.V  = [mpD.V(jp)   ; 0.5.*mpD.V(ip);0.5.*mpD.V(ip)];
    mpD.m  = [mpD.m(jp)   ; 0.5.*mpD.m(ip);0.5.*mpD.m(ip)];
    
    mpD.w  = [ mpD.w(jp,:);[ mpD.w(ip,:); mpD.w(ip,:)]];
    mpD.dD = [mpD.dD(jp,:);[mpD.dD(ip,:);mpD.dD(ip,:)]];
    mpD.de = [mpD.de(:,jp),[mpD.de(:,ip),mpD.de(:,ip)]];
    mpD.ep = [mpD.ep(:,jp),[mpD.ep(:,ip),mpD.ep(:,ip)]];
    mpD.s  = [ mpD.s(:,jp),[ mpD.s(:,ip), mpD.s(:,ip)]];
    mpD.sn = [mpD.sn(:,jp),[mpD.sn(:,ip),mpD.sn(:,ip)]];
    phin= repmat(mpD.phi(ip),2,1);
    mpD.phi = [mpD.phi(jp) phin(:)' ];
    cpn = repmat(mpD.c(ip),2,1);
    mpD.c  = [mpD.c(jp) cpn(:)' ];
    
    mpD.epV  = [ mpD.epV(jp,1);[ mpD.epV(ip,1); mpD.epV(ip,1)]];
    mpD.devep = [ mpD.devep(:,jp),[ mpD.devep(:,ip), mpD.devep(:,ip)]];
    mpD.epII  = [ mpD.epII(1,jp) [ mpD.epII(1,ip)  mpD.epII(1,ip)]];
end

ip=find(mpD.l(:,2)> l.*mpD.l0(:,2));
jp=find(mpD.l(:,2)<=l.*mpD.l0(:,2));
np=2*length(ip)+length(jp);

if(isempty(ip)<1)
    lpx = mpD.l(ip,1).*[1 1];
    lpy = 0.5.*mpD.l(ip,2).*[-1 1];
    
    mpD.x(1:np,1)  = [mpD.x(jp,1)   ;reshape([mpD.x(ip,1) mpD.x(ip,1)],2*size(ip,1),1)     ];
    mpD.x(1:np,2)  = [mpD.x(jp,2)   ;reshape(mpD.x(ip,2)+lpy,2*size(ip,1),1)     ];
    mpD.v(1:np,1)  = [mpD.v(jp,1)   ; reshape(repmat(mpD.v(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.v(1:np,2)  = [mpD.v(jp,2)   ; reshape(repmat(mpD.v(ip,2),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,1) = [mpD.p(jp,1) ; reshape(repmat(mpD.p(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,2) = [mpD.p(jp,2) ; reshape(repmat(mpD.p(ip,2),1,2) ,2*size(ip,1),1)];
    
    mpD.n0 = [mpD.n0(jp)  ; reshape(repmat(mpD.n0(ip),1,2),2*size(ip,1),1)];
    mpD.l  = [ mpD.l(jp,:); [abs(lpx(:)) abs(lpy(:))]                  ];
    mpD.l0 = [mpD.l0(jp,:); [abs(lpx(:)) abs(lpy(:))]                  ];
    mpD.V  = [mpD.V(jp)   ; 0.5.*mpD.V(ip);0.5.*mpD.V(ip)];
    mpD.m  = [mpD.m(jp)   ; 0.5.*mpD.m(ip);0.5.*mpD.m(ip)];

    mpD.w  = [ mpD.w(jp,:);[ mpD.w(ip,:); mpD.w(ip,:)]];
    mpD.dD = [mpD.dD(jp,:);[mpD.dD(ip,:);mpD.dD(ip,:)]];
    mpD.de = [mpD.de(:,jp),[mpD.de(:,ip),mpD.de(:,ip)]];
    mpD.ep = [mpD.ep(:,jp),[mpD.ep(:,ip),mpD.ep(:,ip)]];
    mpD.s  = [ mpD.s(:,jp),[ mpD.s(:,ip), mpD.s(:,ip)]];
    mpD.sn = [mpD.sn(:,jp),[mpD.sn(:,ip),mpD.sn(:,ip)]];
    phin= repmat(mpD.phi(ip),2,1);
    mpD.phi = [mpD.phi(jp) phin(:)' ];
    cpn = repmat(mpD.c(ip),2,1);
    mpD.c  = [mpD.c(jp) cpn(:)' ];
    
        mpD.epV  = [ mpD.epV(jp,1);[ mpD.epV(ip,1); mpD.epV(ip,1)]];
    mpD.devep = [ mpD.devep(:,jp),[ mpD.devep(:,ip), mpD.devep(:,ip)]];
    mpD.epII  = [ mpD.epII(1,jp) [ mpD.epII(1,ip)  mpD.epII(1,ip)]];
end


     



end

