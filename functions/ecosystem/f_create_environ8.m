%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
%
% Creates torus
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [environ] = f_create_environ8 (m,n)

% Determine the surrounding cells of a given cell for calculating
% its environment (coral and algal covers).
% Currently implemented for identifying the 4 surrounding cells
% but can be easily extended to 8, ...
% Note: x gives row indices, y column indices
% For a given cell i, the surrounding cells are defined as follow:
%  			  y
%  		 _ _ _ _ _ _
%  		|_|_|_|_|_|_| 
%  		|_|_|t|_|_|_|  
%  	x	|_|l|i|r|_|_| 
%  		|_|_|b|_|_|_| 
%  		|_|_|_|_|_|_|


%  			  y
%  		 _ _ _ _ _ _
%  		|_|_|_|_|_|_| 
%  		|_|2|3|4|_|_|  
%  	x	|_|9|1|5|_|_| 
%  		|_|8|7|6|_|_| 
%  		|_|_|_|_|_|_|
%
% t: top, b: bottom, l:left, r:right

% Pre-allocate matrix of linear indices (cell identifiers)
environ = uint16(zeros(m*n,9)) ;
% first column stores the cell identifiers of the grid 
environ(:,1) = 1:m*n ;
% extract x- and y-coordinates of every cell
[x,y] = ind2sub([m n], environ(:,1)) ;

% This implements the wrapping of the grid into a torus
%  c1 = [x y] ;
 c2 = [x-1 y-1] ;
 c3 = [x-1 y] ;
 c4 = [x-1 y+1] ;
 c5 = [x y+1] ;
 c6 = [x+1 y+1] ;
 c7 = [x+1 y] ;
 c8 = [x+1 y-1] ;
 c9 = [x y-1] ;

 A REVOIR!!
c3(c3==0) = m ;
c4(c4==0) = m ;
c6(c6==0) = n ;
c8(c8==0) = n ;
c9(c9==0) = n ;
c4(c4 > n) = 1 ;
c5(c5 > n) = 1 ;

c6(c6 > n) = 1 ;
c7(c7 > m) = 1 ;
c8(c8 > m) = 1 ;

c2(c2 > n) = 1 ;
c2(c2==0) = n ;


% This converts to linear indices to be stored in 'environ'
environ(:,2) = sub2ind([m n], c3(:,1), c3(:,2)) ; 
environ(:,3) = sub2ind([m n], c4(:,1), c4(:,2)) ; 
environ(:,4) = sub2ind([m n], c5(:,1), c5(:,2)) ;
environ(:,5) = sub2ind([m n], c6(:,1), c6(:,2)) ; 
environ(:,6) = sub2ind([m n], c7(:,1), c7(:,2)) ; 
environ(:,7) = sub2ind([m n], c8(:,1), c8(:,2)) ; 
environ(:,8) = sub2ind([m n], c9(:,1), c9(:,2)) ; 
environ(:,9) = sub2ind([m n], c2(:,1), c2(:,2)) ; 



