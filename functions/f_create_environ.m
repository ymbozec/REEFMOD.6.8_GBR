%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
%
% Creates torus
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [environ] = f_create_environ (m,n)

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
environ = uint16(zeros(m*n,5)) ;
% first column stores the cell identifiers of the grid 
environ(:,1) = 1:m*n ;
% extract x- and y-coordinates of every cell
[x,y] = ind2sub([m n], environ(:,1)) ;

% This implements the wrapping of the grid into a torus
%  c1 = [x y] ;
%  c2 = [x-1 y-1] ;
%  c3 = [x-1 y] ;
%  c4 = [x-1 y+1] ;
%  c5 = [x y+1] ;
%  c6 = [x+1 y+1] ;
%  c7 = [x+1 y] ;
%  c8 = [x+1 y-1] ;
%  c9 = [x y-1] ;


top = [x-1 y] ; % top cells (c3)
bottom = [x+1 y] ; % bottom cells (c7)
left = [x y-1] ; % left cells (c9)
right = [x y+1] ; % right cells (c5)

top(top==0) = m ;
bottom(bottom > m) = 1 ;
left(left==0) = n ;
right(right > n) = 1 ;

% This converts to linear indices to be stored in 'environ'
environ(:,2) = sub2ind([m n], top(:,1), top(:,2)) ; % top cells
environ(:,3) = sub2ind([m n], bottom(:,1), bottom(:,2)) ; % bottom cells
environ(:,4) = sub2ind([m n], left(:,1), left(:,2)) ; % left cells
environ(:,5) = sub2ind([m n], right(:,1), right(:,2)) ; % right cells