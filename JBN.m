function [Swo, krwJBN, kroJBN] = JBN(Sw_aver, PV_inj, dP, qw_prod, dPo, ...
    muw, muo, krw, kro)
% An implementation of the Johnson–Bossler–Naumann (JBN) method for estimating 
% the relative permeabilities from unsteady-state core flooding data.
%
% SYNOPSIS:
%   [Swo, krwJBN, kroJBN] = JBN(Sw_aver, PV_inj, dP, qw_prod, dPo, muw, muo)
% 
% PARAMETERS:
%   Sw_aver    - An array of average water saturations in the sample.
%   PV_inj     - An array of pore volumes of water, injected in the sample.
%   dP         - An array of differential pressures (in bar)
%   qw_prod    - An array of water production rates (only used to identify 
%                the water breakthrough, i.e. units not used).
%   dPo        - The pressure drop, which corresponds to the single-phase 
%                permeability to oil with at irreducible water saturation
%                (in bar).
%   muw, muo   - Water and oil viscosities (in Pa*sec).
% 
% RETURNS:
%   Swo        - An array of water saturations at the outlet.
%   krwJBN     - The JBN water relative permeability at Swo
%   kroJBN     - The JBN oil relative permeability at Swo

%{
Copyright 2024 Nikolai Andrianov, nia@geus.dk
               Wael Fadi Al-Masri, wfa@geus.dk

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%}


% Consider only the data after water breakthrough
ibt = find(qw_prod > eps, 1);
Sw_aver([1:ibt - 1], :) = [];
PV_inj([1:ibt - 1], :) = [];
dP([1:ibt - 1], :) = [];

% fo/kro = 1 / (muo * lt). Note that the numerical differentiation 
% can yield unphysical values of fo/kro 
fokro = (1 / dPo) * diff(dP ./ PV_inj) ./ diff(1 ./ PV_inj);

% Physical boundaries for fo/kro
min_fokro = min(1/kro, muw/muo/krw);
max_fokro = max(1/kro, muw/muo/krw);

% Remove the unphysical values of fo/kro and cut the corresponding entries in 
% other arrays
ind = fokro >= min_fokro & fokro <= max_fokro;
Sw_aver = Sw_aver(ind);
PV_inj = PV_inj(ind);
fokro = fokro(ind);

% Oil fractional flow at the outlet. Note that the numerical differentiation 
% can yield unphysical values of fo outside of [0, 1]
fo = diff(Sw_aver) ./ diff(PV_inj);

% Remove the unphysical values of fo and cut the corresponding entries in 
% other arrays
ind = fo >= 0 & fo <= 1;
fo = fo(ind);
Sw_aver = Sw_aver(ind);
PV_inj = PV_inj(ind);
fokro = fokro(ind);

% Water saturation at the outlet
Swo = Sw_aver - (fo .* PV_inj);

% krw/kro = fw/fo * muw/muo
kratio = (1 - fo) ./ fo * muw/muo;

kroJBN = fo ./ fokro;
krwJBN = kroJBN .* kratio;

end
