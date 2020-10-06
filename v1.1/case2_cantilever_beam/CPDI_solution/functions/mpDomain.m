function[mpD] = mpDomain(l,mpD);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


lx = 0.5.*sqrt(mpD.r1(:,1).^2+mpD.r1(:,2).^2);
ip = find(lx> l.*mpD.l0(:,1));
jp = find(lx<=l.*mpD.l0(:,1));


np=2*length(ip)+length(jp);

if(isempty(ip)<1)
    lpx = 0.5.*lx(ip).*[-1 1];
    lpy = mpD.l0(ip,2).*[1 1];
    
    mpD.r2(1:np,1)  = [mpD.r2(jp,1)   ; reshape(repmat(     mpD.r2(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.r2(1:np,2)  = [mpD.r2(jp,2)   ; reshape(repmat(     mpD.r2(ip,2),1,2) ,2*size(ip,1),1)];
    mpD.r1(1:np,1)  = [mpD.r1(jp,1)   ; reshape(repmat(0.5.*mpD.r1(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.r1(1:np,2)  = [mpD.r1(jp,2)   ; reshape(repmat(0.5.*mpD.r1(ip,2),1,2) ,2*size(ip,1),1)];
    
    n  = sqrt(mpD.r1(1:np,1).^2+mpD.r1(1:np,2).^2);
    nx = mpD.r1(:,1)./n;
    ny = mpD.r1(:,2)./n;
    
    size(nx)
    size(mpD.r1(:,1))
%     mpD.x(1:np,1) = [mpD.x(jp,1); reshape(mpD.x(ip,1)+lpx,2*size(ip,1),1)     ];
%     mpD.x(1:np,2) = [mpD.x(jp,2); reshape([mpD.x(ip,2) mpD.x(ip,2)],2*size(ip,1),1)     ];
%     xp = repmat(mpD.x(ip,1),1,2)+repmat([0.5 -0.5],length(ip),1).*mpD.r1(ip,1);
%     yp = repmat(mpD.x(ip,2),1,2)+repmat([0.5 -0.5],length(ip),1).*mpD.r1(ip,2);

    xp = repmat(mpD.x(ip,1),1,2)+repmat([-0.25 0.25],length(ip),1).*(nx(ip,1).*n(ip));
    yp = repmat(mpD.x(ip,2),1,2)+repmat([-0.25 0.25],length(ip),1).*(ny(ip,1).*n(ip));
    
    
    
    mpD.x(1:np,1) = [mpD.x(jp,1); reshape(xp,2*size(ip,1),1)];
    mpD.x(1:np,2) = [mpD.x(jp,2); reshape(yp,2*size(ip,1),1)];
    
    mpD.v(1:np,1) = [mpD.v(jp,1); reshape(repmat(mpD.v(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.v(1:np,2) = [mpD.v(jp,2); reshape(repmat(mpD.v(ip,2),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,1) = [mpD.p(jp,1); reshape(repmat(mpD.p(ip,1),1,2) ,2*size(ip,1),1)];
    mpD.p(1:np,2) = [mpD.p(jp,2); reshape(repmat(mpD.p(ip,2),1,2) ,2*size(ip,1),1)];
    
    mpD.n0 = [mpD.n0(jp)  ; reshape(repmat(mpD.n0(ip),1,2),2*size(ip,1),1)];
    mpD.l0 = [mpD.l0(jp,:); [abs(lpx(:)) abs(lpy(:))]                  ];
    mpD.V  = [mpD.V(jp)   ; 0.5.*mpD.V(ip);0.5.*mpD.V(ip)];
    mpD.m  = [mpD.m(jp)   ; 0.5.*mpD.m(ip);0.5.*mpD.m(ip)];
    
    mpD.w  = [ mpD.w(jp,:);[ mpD.w(ip,:); mpD.w(ip,:)]];
    mpD.dD = [mpD.dD(jp,:);[mpD.dD(ip,:);mpD.dD(ip,:)]];
    mpD.F  = [mpD.F(jp,:) ;[mpD.F(ip,:) ;mpD.F(ip,:) ]];
    mpD.e = [mpD.e(:,jp),[mpD.e(:,ip),mpD.e(:,ip)]];
    mpD.ep = [mpD.ep(:,jp),[mpD.ep(:,ip),mpD.ep(:,ip)]];
    mpD.s  = [ mpD.s(:,jp),[ mpD.s(:,ip), mpD.s(:,ip)]];
    mpD.sn = [mpD.sn(:,jp),[mpD.sn(:,ip),mpD.sn(:,ip)]];
    phin= repmat(mpD.phi(ip),2,1);
    mpD.phi = [mpD.phi(jp) phin(:)' ];
    cpn = repmat(mpD.coh(ip),2,1);
    mpD.coh  = [mpD.coh(jp) cpn(:)' ];
    
    mpD.epV  = [ mpD.epV(jp,1);[ mpD.epV(ip,1); mpD.epV(ip,1)]];
    mpD.devep = [ mpD.devep(:,jp),[ mpD.devep(:,ip), mpD.devep(:,ip)]];
    mpD.epII  = [ mpD.epII(1,jp) [ mpD.epII(1,ip)  mpD.epII(1,ip)]];



end

% % 
% % ly = 0.0005.*sqrt(mpD.r2(:,1).^2+mpD.r2(:,2).^2);
% % ip = find(ly> l.*mpD.l0(:,2));
% % jp = find(ly<=l.*mpD.l0(:,2));
% % 
% % 
% % np=2*length(ip)+length(jp);
% % 
% % if(isempty(ip)<1)
% %     lpx = 2.*mpD.l0(ip,1).*[1 1];
% %     lpy = 0.5.*ly(ip).*[-1 1];
% %     
% %     mpD.r2(1:np,1)  = [mpD.r2(jp,1)   ; reshape(repmat(0.5.*mpD.r2(ip,1),1,2) ,2*size(ip,1),1)];
% %     mpD.r2(1:np,2)  = [mpD.r2(jp,2)   ; reshape(repmat(0.5.*mpD.r2(ip,2),1,2) ,2*size(ip,1),1)];
% %     mpD.r1(1:np,1)  = [mpD.r1(jp,1)   ; reshape(repmat(     mpD.r1(ip,1),1,2) ,2*size(ip,1),1)];
% %     mpD.r1(1:np,2)  = [mpD.r1(jp,2)   ; reshape(repmat(     mpD.r1(ip,2),1,2) ,2*size(ip,1),1)];
% % %     mpD.x(1:np,1)  = [mpD.x(jp,1)   ;reshape([mpD.x(ip,1) mpD.x(ip,1)],2*size(ip,1),1)     ];
% % %     mpD.x(1:np,2)  = [mpD.x(jp,2)   ;reshape(mpD.x(ip,2)+lpy,2*size(ip,1),1)     ];
% %     
% % %     xp = (mpD.x(ip,1)+[0.5 -0.5].*mpD.r2(ip,1));
% % %     yp = (mpD.x(ip,2)+[0.5 -0.5].*mpD.r2(ip,2));
% % %     mpD.x(1:np,1) = [mpD.x(jp,1); xp(:)];
% % %     mpD.x(1:np,2) = [mpD.x(jp,2); yp(:)];
% %     xp = repmat(mpD.x(ip,1),1,2)+repmat([-0.5 0.5],length(ip),1).*mpD.r2(ip,1);
% %     yp = repmat(mpD.x(ip,2),1,2)+repmat([-0.5 0.5],length(ip),1).*mpD.r2(ip,2);
% %     mpD.x(1:np,1) = [mpD.x(jp,1); reshape(xp,2*size(ip,1),1)];
% %     mpD.x(1:np,2) = [mpD.x(jp,2); reshape(yp,2*size(ip,1),1)];
% %     
% %     mpD.v(1:np,1)  = [mpD.v(jp,1)   ; reshape(repmat(mpD.v(ip,1),1,2) ,2*size(ip,1),1)];
% %     mpD.v(1:np,2)  = [mpD.v(jp,2)   ; reshape(repmat(mpD.v(ip,2),1,2) ,2*size(ip,1),1)];
% %     mpD.p(1:np,1) = [mpD.p(jp,1) ; reshape(repmat(mpD.p(ip,1),1,2) ,2*size(ip,1),1)];
% %     mpD.p(1:np,2) = [mpD.p(jp,2) ; reshape(repmat(mpD.p(ip,2),1,2) ,2*size(ip,1),1)];
% %     
% %     mpD.n0 = [mpD.n0(jp)  ; reshape(repmat(mpD.n0(ip),1,2),2*size(ip,1),1)];
% %     mpD.l0 = [mpD.l0(jp,:); [abs(lpx(:)) abs(lpy(:))]                  ];
% %     mpD.V  = [mpD.V(jp)   ; 0.5.*mpD.V(ip);0.5.*mpD.V(ip)];
% %     mpD.m  = [mpD.m(jp)   ; 0.5.*mpD.m(ip);0.5.*mpD.m(ip)];
% % 
% %     mpD.w  = [ mpD.w(jp,:);[ mpD.w(ip,:); mpD.w(ip,:)]];
% %     mpD.dD = [mpD.dD(jp,:);[mpD.dD(ip,:);mpD.dD(ip,:)]];
% %     mpD.F  = [mpD.F(jp,:) ;[mpD.F(ip,:) ;mpD.F(ip,:) ]];
% %     mpD.e = [mpD.e(:,jp),[mpD.e(:,ip),mpD.e(:,ip)]];
% %     mpD.ep = [mpD.ep(:,jp),[mpD.ep(:,ip),mpD.ep(:,ip)]];
% %     mpD.s  = [ mpD.s(:,jp),[ mpD.s(:,ip), mpD.s(:,ip)]];
% %     mpD.sn = [mpD.sn(:,jp),[mpD.sn(:,ip),mpD.sn(:,ip)]];
% %     phin= repmat(mpD.phi(ip),2,1);
% %     mpD.phi = [mpD.phi(jp) phin(:)' ];
% %     cpn = repmat(mpD.coh(ip),2,1);
% %     mpD.coh  = [mpD.coh(jp) cpn(:)' ];
% %     
% %         mpD.epV  = [ mpD.epV(jp,1);[ mpD.epV(ip,1); mpD.epV(ip,1)]];
% %     mpD.devep = [ mpD.devep(:,jp),[ mpD.devep(:,ip), mpD.devep(:,ip)]];
% %     mpD.epII  = [ mpD.epII(1,jp) [ mpD.epII(1,ip)  mpD.epII(1,ip)]];
% %     
% % 
% %     
% % end
    
mpD.xc = repmat(mpD.x(:,1),1,4)+mpD.r1(:,1).*[-0.5 0.5 0.5 -0.5]+...
    mpD.r2(:,1).*[-0.5 -0.5 0.5 0.5];
mpD.yc = repmat(mpD.x(:,2),1,4)+mpD.r1(:,2).*[-0.5 0.5 0.5 -0.5]+...
    mpD.r2(:,2).*[-0.5 -0.5 0.5 0.5];
end

