%% Custom ZEMAX Lens Surfaces from MATLAB Data
%
% Replicating "General formula for bi-aspheric singlet lens design free of
% spherical aberration" by R. G. González-Acuña & H. A. Chaparro-Romo
%
% H. Azzouz
% 2021-10-23
% haz@mit.edu

%% Define system configurations
% Settings below are the same as Appendix A
s1 =   1;
n  =   1.5;
t  =   8;
ta = -60;
tb =  70;
basez = ta;
topez = tb + t;
rmax  = 16.8; % Iterative value

ra = linspace(-rmax, rmax, 513); %Original ZEMAX example used 513 points

%% Define first surface function
za = cos(ra/2);
%za = besselj(0, ra);

dra = max(diff(ra)); % should = min(diff(ra))
dza = [diff(za), za(1)-za(end)] ./ dra;

%% Define variables of Eq. (5)
Phi= sqrt(1 - ...
          (ra + ( - ta + za) .* dza) .^ 2 ...
          ./ (n .^ 2 .* (ra .^ 2 + (ta - za) .^ 2) .* (1 + dza .^ 2) ) ...
         ) ...
     ./ sqrt(1 + dza .^ 2);
ri = ( ra + (- ta + za) .* dza ) ...
     ./ ( n .* sqrt(ra .^ 2 + (ta - za) .^ 2) .* (1 + dza .^ 2) ) - ...
     dza .* Phi;
zi = ( dza .* (ra + (- ta + za) .* dza)) ...
     ./ ( n .* sqrt(ra .^ 2 + (ta - za) .^ 2) .* (1 + dza .^ 2) ) + ...
     Phi;
fi = ta - ...
     tb - ...
     sign(ta) .* sqrt( ra .^ 2 + (ta - za) .^ 2 );
h0 = ri .^ 2 .* za + ...
     fi .* n .* zi - ...
     ra .* ri .* zi + ...
     (t + tb) .* zi .^ 2 - ...
     n .^ 2 .* (za + t .* zi);
h1 = ra .^ 2 + ...
     2 .* ra .* ri .* t + ...
     (tb - za) .^ 2 + ...
     t .^ 2 .* (ri .^ 2 + (- 1 + zi) .^ 2) - ...
     2 .* t .* (tb - za) .* (- 1 + zi);
zb = (h0 + ...
      s1 .* sqrt(zi .^ 2 .* (fi .^ 2 - ...
                             2 .* fi .* n .* (ra .* ri + ...
                                              ri .^ 2 .* t + ...
                                              zi .* (t .* (zi - 1) - ...
                                                     tb + za) ...
                                             ) + ...
                             h1 .* n .^ 2 - ...
                             (ra .* zi + ri .* (t + tb - za)) .^ 2 ...
                            ) ...
                ) ...
     ) ...
     ./ (- n .^ 2 + 1);
rb = ra + ...
     (ri .* (- za + zb)) ./ zi;

%% Put B surface on same grid as A surface
zb_new = interpn(rb,zb,ra);
%% System plot
fill([ra fliplr(rb)],[za fliplr(zb)],3/4*[1 1 1])
hold on
plot(ra,za,'k.-')
plot(rb,zb,'b.-')
plot(ra,zb_new,'r.-')
hold off
axis equal tight
xlabel('r  [mm]')
ylabel('z  [m]]')
grid on
title('Replicating Lens From Fig. 4. (d)')
legend('','Surface A','Surface B original',...
       'Surface B regridded','location','best')
%whitebg('black'), set(gcf, 'InvertHardCopy', 'off'); %whitebg('white')

%% Adapted from ZEMAX "How to write a Grid Sag DAT file programmatically"
% https://support.zemax.com/hc/en-us/articles/1500005578182
%
% "In this file, we will use the parameters of surface 1 within the
% downloaded singlet file attached to the KBA to generate a DAT file. The 
% DAT file will produce a sag profile based on the standard sag equation."

%% Intended resolution (npix) and point spacing (dx)
npix = length(ra);
dx = max(diff(ra)); % should = min(diff(ra))

%% Use the arrays za or zb to generate points to place within sag matrix
sagMatrixa = repmat(za,npix,1);
sagMatrixb = repmat(zb,npix,1);

%% Check the sag profile
figure;
surf(sagMatrixa), hold on
surf(sagMatrixb), hold off
xlabel('x'), ylabel('y')
colormap turbo, shading interp
axis square;

%% Flatten the sag matrix
SurfDatClean_a = sagMatrixa(:);
SurfDatClean_b = sagMatrixb(:);

%% Find the Zemax {DATA} folder and create a path to the
% {Zemax}\Objects\Grid Files folder
zemaxData = winqueryreg('HKEY_CURRENT_USER', 'Software\Zemax', 'ZemaxRoot');
gridDir = '\Objects\Grid Files\';

%% Generate the DAT file  -  should functionalize
filename_a = fullfile(zemaxData, gridDir, 'sagData_a.dat');
filename_b = fullfile(zemaxData, gridDir, 'sagData_b.dat');

fileID = fopen(filename_a, 'wt');
fprintf(fileID, '%d %d %1.10g %1.10g %1.f %1.f %1.f', npix, npix, dx, dx, 0, 0, 0);
fprintf(fileID, '\n');
for i = 1:size(SurfDatClean_a)%-1)
    fprintf(fileID, '%0.20g %d %d %d %d\n', SurfDatClean_a(i), 0, 0, 0, 0);
end
fclose(fileID);

fileID = fopen(filename_b, 'wt');
fprintf(fileID, '%d %d %1.10g %1.10g %1.f %1.f %1.f', npix, npix, dx, dx, 0, 0, 0);
fprintf(fileID, '\n');
for i = 1:size(SurfDatClean_b)%-1)
    fprintf(fileID, '%0.20g %d %d %d %d\n', SurfDatClean_b(i), 0, 0, 0, 0);
end
fclose(fileID);

%% end