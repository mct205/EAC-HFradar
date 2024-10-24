function spl = spectrelisse2(x,f,R,p)
% Calcule le spectre brut d'une serie temporelle 
%
% F. G. Schmitt, mai 2006
% O. DURAN MEDINA, octobre 2014
%
% entrée:  - serie x (vecteur)
%          - f = frequence d'echantillonnage
%          - n = taille fenetre FFT
%          - R = nombre de points par decade pour le lissage
%
% sortie:  - un spectre sous la forme (f,E(f))

% Do linear interpolate gaps
indnan = find(isnan(x));
indnotnan = find(~isnan(x));
x(indnan) = interp1(indnotnan,x(indnotnan),indnan);

%nx = max(size(x));
nx = size(x,1);
na = 16;
w = hanning(floor(nx/na));

[spectra,frequency,spectraxxc]= periodogram(x,hamming(nx),nx,f,'ConfidenceLevel',p);

% puis on ne conserve que les premiers points, de 1 a np/2+1

imax = length(frequency);
jmax = fix(R*log10(imax+1)); %arrondi du résultat à l'entier inférieur

pp = zeros(1,jmax); 
ppxxc = zeros(2,jmax); 
nn = zeros(1,jmax); %vecteurs ligne
ff = zeros(1,jmax);

% lissage 
for i = 1:imax ;
    j = fix(R*log10(i+1)); %frequence initiale, creation d'une boite
    pp(j) = spectra(i) + pp(j); %somme des spectres
    ppxxc(1,j) = spectraxxc(i,1) + ppxxc(1,j);
    ppxxc(2,j) = spectraxxc(i,2) + ppxxc(2,j);
    nn(j) = nn(j)+1; %on compte le nombre de frequences dans la boite
    ff(j) = frequency(i)+ff(j); %somme des frequences (moyenne)
end

%Normalisation 
jr = 0; %initialisation
for jj = 1:jmax ;
    
    if (nn(jj) ~=0)
        jr = jr+1;
        spl(jr,2) = pp(jj)/nn(jj);
        spl(jr,3) = ppxxc(1,jj)/nn(jj);
        spl(jr,4) = ppxxc(2,jj)/nn(jj);
        spl(jr,1) = ff(jj)/nn(jj);
        
    end;
   
end

end

