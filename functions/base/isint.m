% ISINT  ( MatLinks) Tests whether or not a number is an integer.
%
%    ISINT(X) returns 1 (true) if X is an integer, or 0 (false) if it is not.
%    If X is a vector, the test will be performed on each element of X and
%    returned as a vector of ones and zeroes corresponding to the results of
%    the test on each respective element of X.
%
%    See also ISEVEN, ISODD, ISNUM, ISNUMERIC.
%
%    Type HELP MATLINKS for a full listing of all MatLinks ToolChest functions.
%
function bool = isint(x)
%===============================================================================
%  Copyright  1999,2000 Julian Andrew de Marchi, Ph.D. (julian@matlinks.net)
%  Use & distribution covered by GNU General Public License (www.gnu.org)
%===============================================================================

bool = (x == fix(x));


%===============================================================================
%  End-Of-File
%===============================================================================
