%sofi=sofiAllFlattenMC(sofi,orders,stddev)
%---------------------------------------
%
%Flatten raw cumulants.
%
%Inputs:
% sofi      Cell array of raw cumulants {Channel 1; Virtual chanel; Channel 2}
% order     Cumulant order of the cumulants
% stddev    Flatten standard deviations if set  {0:flatten mean}
%
%Output:
% sofi      Flat cumulants

%Copyright © 2012 Marcel Leutenegger et al, École Polytechnique Fédérale de Lausanne,
%Laboratoire d'Optique Biomédicale, BM 5.142, Station 17, 1015 Lausanne, Switzerland.

%Adapted for Multi-Color SOFI by Adrien Descloux 29.10.2018

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function sofi=sofiAllFlattenMC(sofi,order,stddev)

if nargin < 3 || isempty(stddev) || ~stddev
    stddev=@mean;
else
    stddev=@std;
end
for channel = 1:length(sofi)
    image=(sofi{channel});
    [X,Y,Nw]=size(image);
    term=ones(order);
    for x=1:order
        for y=1:order
            t=image(x:order:X,y:order:Y,:);
            term(x,y)=feval(stddev,t(:));
        end
    end
    term=mean(term(:))./term;
    term=repmat(term,ceil([X Y]/order));
    sofi{channel}=image(:,:,:).*term(1:X,1:Y,ones(1,Nw));
end
% eof