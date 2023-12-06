% generate toy 2D molecules and their motions to identify functional modes
% 
% INPUT
% Internal creation of molecule and its mutants.
% Three regions:  Outer Left     Inner Middle    Outer Right
% Four geometries:  Free, Linear, Triagonal, Square, Extended 
% Note: Extended is a linear geometry, but in an extended fashion.
% Note: F => Free, L => Linear, T => Triagonal, S => Square, E => Extended
% Assignments of color:
% F => black,  L => Red,  T => Green,  S => Magenta,  E => cyan
% 
% Labeling scheme:
% FFF  FLF  FTF  FSF    EFF  ELF  ETF  ESF   
% FFL  FLL  FTL  FSL    EFL  ELL  ETL  ESL
% FFT  FLT  FTT  FST    EFT  ELT  ETT  EST
%
% Note: FLT => outer left is F,   Inner Middle is L,   outer right is T.
% Note: 2x(3x4) = 12 + 12 = 24 cases to compare. 
% 
% PROCESS
% Simulate the motions using Monte Carlo methods.
% Free => no constraints applied
% Configurations with a local geometry implies constraints are applied. 
% superimpose structures using a home grown method that works well
%% 
clear all             % Start with a clean slate to remove stray variables
%%                                                              User input
nSamples = input('      Enter number of samples: '); 
minNsamples = 4*60;       % 4 samples and must be > 2*na = 58 for 1 sample 
maxNsamples = 999999;
   if( nSamples < minNsamples )
   error(['ERROR: must have at least ',num2str(minNsamples),' samples']);
   end
   if( nSamples > maxNsamples )
   error(['ERROR: must have less than ',num2str(maxNsamples),' samples']);
   end
answ = input(' Enter (y/n) for verification: ','s');
   if( strcmp(answ,'y') )
   reviewCases = true;                                      % review movie
   else
   reviewCases = false;                              % do not review movie
   end
seed = input('Enter seed for random # generator: ');
   if( seed < 1 )
   error('set seed > 1 ');
   end
fname0 = ['A',num2str(nSamples,'%6.6d'),'s',num2str(seed)];
%disp(fname0)
outputFile_eps = 'epsc';                  % eps gives only black and white
outputFile_png = 'png'; 
% error('stop here now');
%%                                                   Initialize parameters
na = 29;
b = 1.2;                                     % bonding length in Angstroms
a = b/1.5;                                   % bumping length in Angstroms
kcb = 500;    % scaled energy/(RT) to get dx ~ 0.05 tolerance in Angstroms
kvdw = 4000;                     % for repulsive van der Wall interactions
nnMax = 4;
nMCshakes = 1000;                           % number of MC steps per shake
nMCsteps0 = 50;                         % to equilibrate, number of shakes
dr0 = 0.01;       % randomize atom positions by 0.01 Angstrom per MC shake
rng(seed);
%%                                                     Initialize molecule
xy0 = zeros(na,2);
nnLink = zeros(na,nnMax);
nnMult = zeros(na,1);
bondList = [ 1  6;  6  7;  7  8;  2  8;  8  9;  9 10;  3 10; 10 11; ...
            11 12;  4 12; 12 13; 13 14;  5 14; 14 15;  6 15; 14 16; ...
            15 16; 16 17; 17 18; 18 19;  7 20; 20 21; 21 22;  9 23; ...
            10 23; 23 24; 24 25; 25 26; 12 27; 27 28; 28 29];
[nb,~] = size(bondList);
   for k=1:nb
   i = bondList(k,1);
   j = bondList(k,2);
   nnMult(i) = nnMult(i) + 1;
   nnMult(j) = nnMult(j) + 1;
      if( nnMult(i) > nnMax )
      error(['Degree of atom ',num2str(i),' is greater than nnMax = ', ...
              num2str(nnMax)]);
      end
      if( nnMult(j) > nnMax )
      error(['Degree of atom ',num2str(j),' is greater than nnMax = ', ...
              num2str(nnMax)]);
      end
   nnLink(i, nnMult(i) ) = j;
   nnLink(j, nnMult(j) ) = i;
   end
 % -------------------------------------------------------- position atoms 
 % From geometry, inscribing a regular decahedron in a circle of radius g.
 % For 10 sides, each wedge has (2*pi)/10 radians or 36 degrees. 
 % As such the formula for g is: g = (b/2)/sin(18 degrees);
 g = (b/2)/sin(pi/10);
 circle1Atoms = [7, 8, 9, 10, 11, 12, 13, 14, 15, 6];
    for i=1:10
    j = circle1Atoms(i);
    phi = 0.5*pi - (pi/5)*(i-1);
    xy0(j,1) = g*cos(phi);                                   % x-component
    xy0(j,2) = g*sin(phi);                                   % y-component
    end
%  figure(1);
%  clf;
%  plot(xy0(circle1Atoms,1),xy0(circle1Atoms,2),'o');
 % -------------------------------------------------------- put inner loop
 circle2Atoms = [2, 3, 4, 5, 1];
 circle2AtomsSign = [-1, 1, 1, -1, 1];
   for i=1:5
   j = circle2Atoms(i);
   s = circle2AtomsSign(i);
   phi = 0.5*pi - (pi/5) - (2*pi/5)*(i-1); 
   xy0(j,1) = (g-s*b)*cos(phi);                              % x-component
   xy0(j,2) = (g-s*b)*sin(phi);                              % y-component
   end
% hold on;
% plot(xy0(circle2Atoms,1),xy0(circle2Atoms,2),'ro');
% pbaspect([1 1 1]);
% -------------------------------------------------- put on vertical wings
xy0(20,1) = 0;
xy0(21,1) = 0;
xy0(22,1) = 0;
xy0(20,2) = g + b;
xy0(21,2) = g + 2*b;
xy0(22,2) = g + 3*b;
% ----------------------
xy0(27,1) = 0;
xy0(28,1) = 0;
xy0(29,1) = 0;
xy0(27,2) = -g - b;
xy0(28,2) = -g - 2*b;
xy0(29,2) = -g - 3*b;
% ------------------------------------------------ put on horizontal wings
xh = g*sin(2*pi/5) + b*sind(60);
xy0(23,1) = xh;
xy0(24,1) = xh + b;
xy0(25,1) = xh + 2*b;
xy0(26,1) = xh + 3*b;
xy0(23,2) = 0;
xy0(24,2) = 0;
xy0(25,2) = 0;
xy0(26,2) = 0;
% ----------------------
xh = -g*sin(2*pi/5) - b*sind(60);
xy0(16,1) = xh;
xy0(17,1) = xh - b;
xy0(18,1) = xh - 2*b;
xy0(19,1) = xh - 3*b;
xy0(16,2) = 0;
xy0(17,2) = 0;
xy0(18,2) = 0;
xy0(19,2) = 0;
%%                                                          create mutants
% -------------------------------
% color code: 0   1   2   3   4
%   geometry: F   L   T   S   E
%      color: k   r   g   m   c
% -------------------------------
nMutantMolecules = 24;
moleculeCases = 1:nMutantMolecules;      % list the cases to be considered
nAddedBonds = cell(1,nMutantMolecules);
addedBonds = cell(1,nMutantMolecules);
bondLengths = cell(1,nMutantMolecules);
mutations = cell(1,nMutantMolecules);
nMutations = cell(1,nMutantMolecules);
nameMolecule = cell(1,nMutantMolecules);
Mname = cell(1,nMutantMolecules);
dataMatrix = cell(1,nMutantMolecules);
Aname = cell(1,nMutantMolecules);
% -------------------------------------------------------------- case: FFF
m = 1; 
nMutations{m} = 0;
mutations{m} = [];
colourize{m} = [];
addedBonds{m} = [];
nAddedBonds{m} = 0;
bondLengths{m} = [];
nameMolecule{m} = 'molecule: FFF';
Mname{m} = 'moleculeFFF';
Aname{m} = [fname0,'FFF'];
% -------------------------------------------------------------- case: FFL
m = 2;
nMutations{m} = 3;
mutations{m} = [2, 22, 26];
colourize{m} = [1,  1,  1];                               % L => red color
addedBonds{m} = [22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FFL';
Mname{m} = 'moleculeFFL';
Aname{m} = [fname0,'FFL'];
% -------------------------------------------------------------- case: FFT
m = 3;
nMutations{m} = 3;
mutations{m} = [2, 22, 26];
colourize{m} = [2,  2,  2];                             % T => green color
addedBonds{m} = [22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FFT';
Mname{m} = 'moleculeFFT';
Aname{m} = [fname0,'FFT'];
% -------------------------------------------------------------- case: FLF
m = 4;
nMutations{m} = 3;
mutations{m} = [1, 3, 4];
colourize{m} = [1, 1, 1];                                 % L => red color
addedBonds{m} = [1,4; 1,3; 4,3];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.5*b;
bx3 = b14/2;
bondLengths{m} = [b14; bx3; bx3];
nameMolecule{m} = 'molecule: FLF';
Mname{m} = 'moleculeFLF';
Aname{m} = [fname0,'FLF'];
% -------------------------------------------------------------- case: FLL
m = 5;
nMutations{m} = 6;
mutations{m} = [1, 3, 4, 2, 22, 26];
colourize{m} = [1, 1, 1, 1,  1,  1];                            % 1 => red
addedBonds{m} = [1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.5*b;
bx3 = b14/2;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FLL';
Mname{m} = 'moleculeFLL';
Aname{m} = [fname0,'FLL'];
% -------------------------------------------------------------- case: FLT
m = 6;
nMutations{m} = 6;
mutations{m} = [1, 3, 4, 2, 22, 26];
colourize{m} = [1, 1, 1, 2,  2,  2];                % (1,2) => (red,green)
addedBonds{m} = [1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.5*b;
bx3 = b14/2;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FLT';
Mname{m} = 'moleculeFLT';
Aname{m} = [fname0,'FLT'];
% -------------------------------------------------------------- case: FTF
m = 7;
nMutations{m} = 3;
mutations{m} = [1, 3, 4];
colourize{m} = [2, 2, 2];                               % T => green color
addedBonds{m} = [1,4; 1,3; 4,3];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.0*b;
bx3 = b14;
bondLengths{m} = [b14; bx3; bx3];
nameMolecule{m} = 'molecule: FTF';
Mname{m} = 'moleculeFTF';
Aname{m} = [fname0,'FTF'];
% -------------------------------------------------------------- case: FTL
m = 8;
nMutations{m} = 6;
mutations{m} = [1, 3, 4, 2, 22, 26];
colourize{m} = [2, 2, 2, 1,  1,  1];                % (2,1) => (green,red)
addedBonds{m} = [1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.0*b;
bx3 = b14;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FTL';
Mname{m} = 'moleculeFTL';
Aname{m} = [fname0,'FTL'];
% -------------------------------------------------------------- case: FTT
m = 9;
nMutations{m} = 6;
mutations{m} = [1, 3, 4, 2, 22, 26];
colourize{m} = [2, 2, 2, 2,  2,  2];                          % 2 => green
addedBonds{m} = [1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b14 = 1.0*b;
bx3 = b14;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FTT';
Mname{m} = 'moleculeFTT';
Aname{m} = [fname0,'FTT'];
% -------------------------------------------------------------- case: FSF
m = 10;
nMutations{m} = 4;
mutations{m} = [1, 3, 4, 14];
colourize{m} = [3, 3, 3,  3];                         % 3 => magenta color
addedBonds{m} = [1,3; 3,4; 4,14; 14,1; 14,3; 1,4];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
bondLengths{m} = [bside; bside; bside; bside; bdiag; bdiag];
nameMolecule{m} = 'molecule: FSF';
Mname{m} = 'moleculeFSF';
Aname{m} = [fname0,'FSF'];
% -------------------------------------------------------------- case: FSL
m = 11;
nMutations{m} = 7;
mutations{m} = [1, 3, 4, 14, 2, 22, 26];
colourize{m} = [3, 3, 3,  3, 1,  1,  1];    % (3,1) => (magenta,red) color
addedBonds{m} = [1,3; 3,4; 4,14; 14,1; 14,3; 1,4; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [bside; bside; bside; bside; bdiag; bdiag; ...
                                               b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FSL';
Mname{m} = 'moleculeFSL';
Aname{m} = [fname0,'FSL'];
% -------------------------------------------------------------- case: FST
m = 12;
nMutations{m} = 7;
mutations{m} = [1, 3, 4, 14, 2, 22, 26];
colourize{m} = [3, 3, 3,  3, 2,  2,  2];  % (3,2) => (magenta,green) color
addedBonds{m} = [1,3; 3,4; 4,14; 14,1; 14,3; 1,4; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [bside; bside; bside; bside; bdiag; bdiag; ...
                                               b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: FST';
Mname{m} = 'moleculeFST';
Aname{m} = [fname0,'FST'];
% -------------------------------------------------------------- case: EFF
m = 13; 
nMutations{m} = 3;
mutations{m} = [19, 5, 29];
colourize{m} = [ 4, 4,  4];                         % E => 4 => cyan color
addedBonds{m} = [19,29; 19,5; 5,29];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
bondLengths{m} = [b19_29; bx5; bx5];
nameMolecule{m} = 'molecule: EFF';
Mname{m} = 'moleculeEFF';
Aname{m} = [fname0,'EFF'];
% -------------------------------------------------------------- case: EFL
m = 14; 
nMutations{m} = 6;
mutations{m} = [19, 5, 29, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 1,  1,  1];         % (4,1) => (cyan,red) color
addedBonds{m} = [19,29; 19,5; 5,29; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b19_29; bx5; bx5; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: EFL';
Mname{m} = 'moleculeEFL';
Aname{m} = [fname0,'EFL'];
% -------------------------------------------------------------- case: EFT
m = 15; 
nMutations{m} = 6;
mutations{m} = [19, 5, 29, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 2,  2,  2];       % (4,2) => (cyan,green) color
addedBonds{m} = [19,29; 19,5; 5,29; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b19_29; bx5; bx5; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: EFT';
Mname{m} = 'moleculeEFT';
Aname{m} = [fname0,'EFT'];
% -------------------------------------------------------------- case: ELF
m = 16; 
nMutations{m} = 6;
mutations{m} = [19, 5, 29, 1, 3, 4];
colourize{m} = [ 4, 4,  4, 1, 1, 1];           % (4,1) => (cyan,red) color
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.5*b;
bx3 = b14/2;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3];
nameMolecule{m} = 'molecule: ELF';
Mname{m} = 'moleculeELF';
Aname{m} = [fname0,'ELF'];
% -------------------------------------------------------------- case: ELL
m = 17; 
nMutations{m} = 9;
mutations{m} = [19, 5, 29, 1, 3, 4, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 1, 1, 1, 1,  1,  1];      % (4,1) => (cyan,red)
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.5*b;
bx3 = b14/2;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: ELL';
Mname{m} = 'moleculeELL';
Aname{m} = [fname0,'ELL'];
% -------------------------------------------------------------- case: ELT
m = 18; 
nMutations{m} = 9;
mutations{m} = [19, 5, 29, 1, 3, 4, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 1, 1, 1, 2,  2,  2];  % (4,1,2) => (cyan,red,g)
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.5*b;
bx3 = b14/2;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: ELT';
Mname{m} = 'moleculeELT';
Aname{m} = [fname0,'ELT'];
% -------------------------------------------------------------- case: ETF
m = 19; 
nMutations{m} = 6;
mutations{m} = [19, 5, 29, 1, 3, 4];
colourize{m} = [ 4, 4,  4, 2, 2, 2];         % (4,2) => (cyan,green) color
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.0*b;
bx3 = b14;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3];
nameMolecule{m} = 'molecule: ETF';
Mname{m} = 'moleculeETF';
Aname{m} = [fname0,'ETF'];
% -------------------------------------------------------------- case: ETL
m = 20; 
nMutations{m} = 9;
mutations{m} = [19, 5, 29, 1, 3, 4, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 2, 2, 2, 1,  1,  1];  % (4,2,1) => (cyan,g,red)
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.0*b;
bx3 = b14;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: ETL';
Mname{m} = 'moleculeETL';
Aname{m} = [fname0,'ETL'];
% -------------------------------------------------------------- case: ETT
m = 21; 
nMutations{m} = 9;
mutations{m} = [19, 5, 29, 1, 3, 4, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 2, 2, 2, 2,  2,  2];    % (4,2) => (cyan,green)
addedBonds{m} = [19,29; 19,5; 5,29; 1,4; 1,3; 4,3; 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
b14 = 1.0*b;
bx3 = b14;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b19_29; bx5; bx5; b14; bx3; bx3; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: ETT';
Mname{m} = 'moleculeETT';
Aname{m} = [fname0,'ETT'];
% -------------------------------------------------------------- case: ESF
m = 22;
nMutations{m} = 7;
mutations{m} = [19, 5, 29, 1, 3, 4, 14];
colourize{m} = [ 4, 4,  4, 3, 3, 3,  3];         % (4,3) => (cyan,magenta)
addedBonds{m} = [19,29; 19,5; 5,29; 1,3; 3,4; 4,14; 14,1; 14,3; 1,4];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
bondLengths{m} = [b19_29; bx5; bx5; bside; bside; bside; bside; ...
                                                  bdiag; bdiag];
nameMolecule{m} = 'molecule: ESF';
Mname{m} = 'moleculeESF';
Aname{m} = [fname0,'ESF'];
% -------------------------------------------------------------- case: ESL
m = 23;
nMutations{m} = 10;
mutations{m} = [19, 5, 29, 1, 3, 4, 14, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 3, 3, 3,  3, 1,  1,  1];   % (4,3,1) => (c,m,r)
addedBonds{m} = [19,29; 19,5; 5,29; 1,3; 3,4; 4,14; 14,1; 14,3; 1,4; ...
                 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
b22_26 = 1.5*b;
bx2 = b22_26/2;
bondLengths{m} = [b19_29; bx5; bx5; bside; bside; bside; bside; ...
                               bdiag; bdiag; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: ESL';
Mname{m} = 'moleculeESL';
Aname{m} = [fname0,'ESL'];
% -------------------------------------------------------------- case: EST
m = 24;
nMutations{m} = 10;
mutations{m} = [19, 5, 29, 1, 3, 4, 14, 2, 22, 26];
colourize{m} = [ 4, 4,  4, 3, 3, 3,  3, 2,  2,  2];   % (4,3,2) => (c,m,g)
addedBonds{m} = [19,29; 19,5; 5,29; 1,3; 3,4; 4,14; 14,1; 14,3; 1,4; ...
                 22,26; 22,2; 26,2];
[nb_add,~] = size(addedBonds{m});
nAddedBonds{m} = nb_add;
b19_29 = 7.5*b;
bx5 = b19_29/2;
bside = 1.0*b;
bdiag = sqrt(2)*bside;
b22_26 = 1.0*b;
bx2 = b22_26;
bondLengths{m} = [b19_29; bx5; bx5; bside; bside; bside; bside; ...
                               bdiag; bdiag; b22_26; bx2; bx2];
nameMolecule{m} = 'molecule: EST';
Mname{m} = 'moleculeEST';
Aname{m} = [fname0,'EST'];
% ---------------------------------------------------------- for debugging
   if( reviewCases )
      for m=1:nMutantMolecules
      disp('     ');
      disp(['case ',num2str(m),'      ',nameMolecule{m}]);
      disp(['nAddedBonds = ',num2str( nAddedBonds{m} )]);
      disp('----------------------------------');
      disp( num2str([addedBonds{m},bondLengths{m}]) );
      disp(['nMutations = ',num2str( nMutations{m} )]);
         if( nMutations{m} > 0 ) 
         disp(['mutations = ',num2str(mutations{m})]);
         end
      % show initial configuration
      gcf = drawMolecule(m,xy0,bondList,mutations{m}, ...
                                        colourize{m},nameMolecule{m});
      pause(3)
      disp('    ');
      end
   end
%%                                          set force field parameters FFP
FFP = struct;
FFP.b = b;
FFP.a = a;
FFP.kcb = kcb;
FFP.kvdw = kvdw;
FFP.na = na;
FFP.nb = nb;
%%                                           initialize alignment template
alignAtomList = 6:15;
%alignAtomList = 1:na;
[xRot,yRot,rotA] = getRotTemplate(xy0,alignAtomList); 
                                % ^^^--> xy0 is used as a common reference
%%                                              run simulations per mutant
Nm = length(moleculeCases);
AddedEnergy = struct;
na2 = 2*na;                     % for data matrix: need x and y components
for mmm=1:Nm           
m = moleculeCases(mmm);
Adata = zeros(na2,nSamples);
AddedEnergy.addedBonds = addedBonds{m};
AddedEnergy.nB = nAddedBonds{m};
AddedEnergy.bondLengths = bondLengths{m};
xy = xy0;
%%                                                      equilibrate system
countMove = 0;
countFail = 0;
energy = getEnergy(xy,bondList,AddedEnergy,FFP);
dr = dr0;
   for frames=1:nMCsteps0
      for shakes=1:nMCshakes
      xy_test = xy + dr*randn(na,2);
      energy_test = getEnergy(xy_test,bondList,AddedEnergy,FFP); 
      dE = energy_test - energy;
         if( dE < 0 )
         xy = xy_test;
         energy = energy_test;
         countMove = countMove + 1;
         else
         prob = exp(-dE);            % RT is already built into the factor
            if( rand < prob )
            xy = xy_test;
            energy = energy_test;
            countMove = countMove + 1;
            else
            countFail = countFail + 1;
            end
         end
      end
% ----------------------------------------- place center of mass at origin
   xyave = mean(xy);
   xy = xy - xyave;
% ---------------------------------------------------- align 2d structures
   xy = align2Dstructure(xRot,yRot,rotA,xy,alignAtomList);
   gcf = drawMolecule(m,xy,bondList,mutations{m}, ...
                                    colourize{m},nameMolecule{m});
   frac = countMove/(countMove + countFail);
      if( frac < 0.15 )            % => try for more than 15% success rate 
      dr = dr*0.95;
      elseif( frac > 0.35 )        % => try for less than 35% success rate
      dr = dr*1.05;
      end
   %disp([frac,dr]);            % => To check success rate and adaptive dr
   pause(0.1);
   end
% ---------------- draw out initial starting structures for production run
%h = drawMolecule(m,xy,bondList,mutations{m}, ...
%colourize{m},nameMolecule{m});
%fname = ['molecule',num2str(m)];
%saveas(h,Mname{m},outputFile_png);
%saveas(h,Mname{m},outputFile_eps);
%%                                                          production run
countMove = 0;
countFail = 0;
energy = getEnergy(xy,bondList,AddedEnergy,FFP);
   for frames=1:nSamples
      for shakes=1:nMCshakes
      xy_test = xy + dr*randn(na,2);
      energy_test = getEnergy(xy_test,bondList,AddedEnergy,FFP); 
      dE = energy_test - energy;
         if( dE < 0 )
         xy = xy_test;
         energy = energy_test;
         countMove = countMove + 1;
         else
         prob = exp(-dE);            % RT is already built into the factor
            if( rand < prob )
            xy = xy_test;
            energy = energy_test;
            countMove = countMove + 1;
            else
            countFail = countFail + 1;
            end
         end
      end
% ----------------------------------------- place center of mass at origin
   xyave = mean(xy);
   xy = xy - xyave;
% ---------------------------------------------------- align 2d structures
   xy = align2Dstructure(xRot,yRot,rotA,xy,alignAtomList);
% -------------------------------- record conformation in the Adata matrix
   Adata(:,frames) = xy(1:end);          % packing all x followed by all y
   %gcf = drawMolecule(m,xy,bondList,mutations{m}, ...
   %                                 colourize{m},nameMolecule{m});
   frac = countMove/(countMove + countFail);
      if( frac < 0.15 )            % => try for more than 15% success rate 
      dr = dr*0.95;
      elseif( frac > 0.35 )        % => try for less than 35% success rate
      dr = dr*1.05;
      end
%    disp([frac,dr]);
%    pause(0.1);
   end
gcf = drawMolecule(m,xy,bondList,mutations{m}, ...
                                 colourize{m},nameMolecule{m});
dataMatrix{m} = Adata;
dlmwrite(Aname{m},dataMatrix{m});
end
%%                                                
         
function gcf_out = drawMolecule(nFig,xy,bondList,mutations, ...
                                                 colourize,namePlot)
% define figure
figure(nFig);
cmap = ['r','g','m','c'];
clf;
% ------------------------------------------------------ place bonds first
hold on;
[nb,~] = size(bondList);
   for k=1:nb
   i = bondList(k,1);
   j = bondList(k,2);
   xBond = [xy(i,1),xy(j,1)];
   yBond = [xy(i,2),xy(j,2)];
   plot(xBond,yBond,'k','linewidth',1.0);
   end
% ------------------------------------------------- plot all generic atoms
plot(xy(:,1),xy(:,2),'o','MarkerEdgeColor','k', ...
                         'MarkerFaceColor','k','MarkerSize',8);
% -------------------------------- overwrite atoms if atoms are functional
Nmutants = length(mutations);
   if( Nmutants > 0 )
   hold on;
      for i=1:4
      L = ( i == colourize );
         if( ~isempty(L) )
         c = cmap(i);
         pts = mutations(L);
         plot(xy(pts,1),xy(pts,2),'o','MarkerEdgeColor',c, ...
                    'MarkerFaceColor',c,'MarkerSize',9);                           
         end
      end
   end
% ------------------------------------------------- make graph look pretty
daspect([1 1 1]);
xlim([-7,7]);
ylim([-7,7]);
xlabel('x');
ylabel('y');
title(namePlot);
gcf_out = gcf;
end

function energy = getEnergy0(xy,bondList,FFP)
% b = bonding length in Angstroms
% a = bumping length in Angstroms
% kcb = scaled energy/(RT) to get dx ~ 0.05 tolerance in Angstroms
% kvdw = scaled energy/RT for repulsive van der Wall interactions
%%                                                  setup local parameters
b = FFP.b;
a = FFP.a;
kcb = FFP.kcb;
kvdw = FFP.kvdw;
na = FFP.na;
nb = FFP.nb;
a2 = a*a;
%%                                 calculate total energy of configuration
energy = 0;
%%                                             covalent bond contributions
   for k=1:nb
   i = bondList(k,1);
   j = bondList(k,2);
   xi = xy(i,1);
   yi = xy(i,2);
   xj = xy(j,1);
   yj = xy(j,2);
   d = sqrt( (xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) ); 
   energy = energy + kcb*(b-d)*(b-d);
   end
%%                                             van der Waals contributions
   for i=1:na
      for j=i+1:na
      xi = xy(i,1);
      yi = xy(i,2);
      xj = xy(j,1);
      yj = xy(j,2);
      d2 = (xi - xj)*(xi - xj) + (yi - yj)*(yi - yj);
         if( d2 < a2 )
         d = sqrt(d2);
         x = a - d;
         energy = energy + kvdw*x*x;
         end
      end
   end
end

function energy = getEnergy(xy,bondList,AddedEnergy,FFP)
% b = bonding length in Angstroms
% a = bumping length in Angstroms
% kcb = scaled energy/(RT) to get dx ~ 0.05 tolerance in Angstroms
% kvdw = scaled energy/RT for repulsive van der Wall interactions
%%                                                  setup local parameters
kcb = FFP.kcb;
%%                                 calculate total energy of configuration
energy = getEnergy0(xy,bondList,FFP);                        % base energy
%                                                     add the added energy
addedBonds = AddedEnergy.addedBonds;
nB = AddedEnergy.nB;
b = AddedEnergy.bondLengths;
%                            for each bond there is an energy contribution
   for k=1:nB
   i = addedBonds(k,1);
   j = addedBonds(k,2);
   x1 = xy(i,1);
   y1 = xy(i,2);
   x2 = xy(j,1);
   y2 = xy(j,2);
   d1 = sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) );
   d2 = b(k);
   energy = energy + kcb*(d2 - d1)*(d2 - d1);
   end
end

function [xRot,yRot,rotA] = getRotTemplate(xy,alignAtomList)
x = xy(alignAtomList,1);
y = xy(alignAtomList,2);
xave = mean(x);
yave = mean(y);
x = x - xave;
y = y - yave;
na = length(x);
rotA = -0.15:0.004:0.15;          % rotation angle for reference structure 
nPhi = length(rotA);
xRot = zeros(na,nPhi);
yRot = zeros(na,nPhi);
    for k=1:nPhi
    dPhi = rotA(k);
    c = cos(dPhi);
    s = sin(dPhi);
    rotM = [ c -s; s c];
       for j=1:na
       result = rotM*[x(j); y(j)];
       xRot(j,k) = result(1);
       yRot(j,k) = result(2);
       end
    end
rotA = -rotA;                           % rotation angle for new structure
end

function xy = align2Dstructure(xRot,yRot,rotA,xy,alignAtomList)
x = xy(alignAtomList,1);
y = xy(alignAtomList,2);
xave = mean(x);
yave = mean(y);
x = x - xave;
y = y - yave;
nPhi = length(rotA);
energy = zeros(1,nPhi);
   for k=1:nPhi
   x0 = xRot(:,k);
   y0 = yRot(:,k);
   energy(k) = sum( (x - x0).^2 ) + sum( (y - y0).^2 ); 
   end
[~,k] = min(energy);
dPhi = rotA(k);
c = cos(dPhi);
s = sin(dPhi);
rotM = [ c -s; s c];
[na,~] = size(xy);
   for j=1:na
   result = rotM*[xy(j,1); xy(j,2)];
   xy(j,1) = result(1);
   xy(j,2) = result(2);
   end
end 