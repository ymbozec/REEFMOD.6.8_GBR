% ISEVEN  ( MatLinks) Tests whether an integer is even or not.
%
%    ISEVEN(X) returns 1 (true) if X is an even integer, or 0 (false) if it is
%    not.  If X is not an integer, 0 is returned.  If X is a vector, the test
%    will be performed on each element of X and returned as a vector of ones
%    and zeroes corresponding to the results of the test on each respective
%    element of X.
%
%    See also ISODD, ISINT, ISNUM, ISNUMERIC.
%
%    Type HELP MATLINKS for a full listing of all MatLinks ToolChest functions.
%
function bool = iseven(x)
%===============================================================================
%  Copyright  1999,2000 Julian Andrew de Marchi, Ph.D. (julian@matlinks.net)
%  Use & distribution covered by GNU General Public License (www.gnu.org)
%===============================================================================

bool = isint(x./2);


%===============================================================================
%  End-Of-File
%===============================================================================
