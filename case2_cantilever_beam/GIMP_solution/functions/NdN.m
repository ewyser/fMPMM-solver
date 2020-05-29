function [mpD] = NdN(meD,mpD,p2N)
% f_NDN : function which calculate the basis function between
% material points and nodes through the array p2N, which is a topological
% indexing between a material point and all the nodes it integrates to
%   Detailed explanation goes here
%% COMPUTE (X,Y)-BASIS FUNCTION
D        = (repmat(mpD.x(:,1),1,meD.nNe) - meD.x(p2N))                    ;% distance x from point to node
[Nx,dNx] = N(D,meD.h(1),repmat(mpD.l(:,1),1,meD.nNe))                     ;% - see function      
D        = (repmat(mpD.x(:,2),1,meD.nNe) - meD.y(p2N))                    ;% distance y from point to node
[Ny,dNy] = N(D,meD.h(2),repmat(mpD.l(:,2),1,meD.nNe))                     ;% - see function 
%-------------------------------------------------------------------------%

%% CONVOLUTION OF BASIS FUNCTIONS
mpD.N    =  Nx.* Ny                                                       ;% basis 
mpD.dNx  = dNx.* Ny                                                       ;% x-gradient basis 
mpD.dNy  =  Nx.*dNy                                                       ;% y-gradient basis  
%-------------------------------------------------------------------------%

%% B MATRIX ASSEMBLY
iDx            = 1:meD.DoF:meD.nDF(1)-1                                   ;% x component global node index
iDy            = iDx+1                                                    ;% y component global node index
mpD.B(1,iDx,:) = mpD.dNx'                                                 ;% -
mpD.B(2,iDy,:) = mpD.dNy'                                                 ;% -
mpD.B(4,iDx,:) = mpD.dNy'                                                 ;% -
mpD.B(4,iDy,:) = mpD.dNx'                                                 ;% -
%-------------------------------------------------------------------------%

end
function [N,dN]=N(dX,h,lp)
%% COMPUTE BASIS FUNCTIONS
lp = 2*lp                                                                 ;% length of mp domain
c1 = ( abs(dX)< (  0.5*lp)                        )                       ;% logical array 1
c2 = ((abs(dX)>=(  0.5*lp)) & (abs(dX)<(h-0.5*lp)))                       ;% logical array 2
c3 = ((abs(dX)>=(h-0.5*lp)) & (abs(dX)<(h+0.5*lp)))                       ;% logical array 3
% BASIS FUNCTION
N1 = 1-((4*dX.^2+lp.^2)./(4*h.*lp))                                       ;% basis function according to c1
N2 = 1-(abs(dX)./h)                                                       ;% basis function according to c2
N3 = ((h+0.5*lp-abs(dX)).^2)./(2*h.*lp)                                   ;% basis function according to c3
N  = c1.*N1+c2.*N2+c3.*N3                                                 ;% basis function
% BASIS FUNCTION GRADIENT
dN1= -((8*dX)./(4*h.*lp))                                                 ;% gradient basis function according to c1
dN2= sign(dX).*(-1/h)                                                     ;% gradient basis function according to c2
dN3=-sign(dX).*(h+0.5*lp-abs(dX))./(h*lp)                                 ;% gradient basis function according to c3
dN = c1.*dN1+c2.*dN2+c3.*dN3                                              ;% gradient basis function
%-------------------------------------------------------------------------%

% % %% COMPUTE BASIS FUNCTIONS
% % c1 = (((-h-lp)<dX) & (dX<=(-h+lp)))                                       ;% logical array 1
% % c2 = (((-h+lp)<dX) & (dX<=(  -lp)))                                       ;% logical array 2
% % c3 = ((    -lp<dX) & (dX<=    lp) )                                       ;% logical array 3
% % c4 = ((     lp<dX) & (dX<=( h-lp)))                                       ;% logical array 4
% % c5 = ((( h-lp)<dX) & (dX<=( h+lp)))                                       ;% logical array 5
% % % BASIS FUNCTION
% % N1 = ((h+lp+dX).^2)./(4*h*lp)                                             ;% basis function according to c1
% % N2 = 1+(dX./h)                                                            ;% basis function according to c2
% % N3 = 1-((dX.^2+lp.^2)./(2*h*lp))                                          ;% basis function according to c3
% % N4 = 1-(dX./h)                                                            ;% basis function according to c4
% % N5 = ((h+lp-dX).^2)./(4*h*lp)                                             ;% basis function according to c5
% % N  = c1.*N1+c2.*N2+c3.*N3+c4.*N4+c5.*N5                                   ;% basis function
% % % BASIS FUNCTION GRADIENT
% % dN1 = (h+lp+dX)./(2*h*lp)                                                 ;% gradient basis function according to c1
% % dN2 = 1/h.*ones(size(dX))                                                 ;% gradient basis function according to c2
% % dN3 = -dX./(h*lp)                                                         ;% gradient basis function according to c3
% % dN4 = -1/h.*ones(size(dX))                                                ;% gradient basis function according to c4
% % dN5 = -(h+lp-dX)./(2*h*lp)                                                ;% gradient basis function according to c5
% % dN  = c1.*dN1+c2.*dN2+c3.*dN3+c4.*dN4+c5.*dN5                             ;% gradient basis function
% % %-------------------------------------------------------------------------%

end