
clear   % clear workspace
clc     % clear console screen

%--------------------------------------------------------------------------

%*************************%
%* Adjustable Parameters *%
%*************************%

% The location of the folder that contains the data
path='dataset1\';

% the different datasets
datasets={'hmdad' 'disbiome'};

% CLASSIFIER -------------------------------
classifier='iPALMGLMF';
% ------------------------------------------

% Parameters and Options -------------------
%WKNKN
use_WKNKN = 1;      % 1=yes, 0=no
K = 5;              % number of K nearest known neighbors
eta = 0.7;          % decay rate (also used by WNN in RLS-WNN)

% CROSS VALIDATION PARAMETERS --------------
n = 5;  % the 'n' in "n-fold experiment"
% m = 450; % number of known associations in the hmdad database
m = 4351; % number of known associations in the disbiome database

% Terminology:
% Y = Interaction matrix
% GM = microbe similarity matrix  // microbe
% GD = disease similarity matrix  //disease

disp('==============================================================');
fprintf('\nClassifier Used: %s',classifier);
fprintf('\nCV Setting Used: 5-Fold\n');

if use_WKNKN
    fprintf('\nusing WKNKN: K=%i, eta=%g\n',K,eta);
end
fprintf('\n');

% Select ds based on the value of m
if m == 450
    ds = 1;
elseif m == 4351
    ds = 2;
else
    error('Unknown value of n');
end
% Process the selected dataset
disp('--------------------------------------------------------------');
fprintf('\nData Set: %s\n', datasets{ds});

% Load data
[Y, GM, GD] = getdata(path, datasets{ds});

% Perform prediction and print evaluation metrics
[y3, matrixs] = crossValidation(Y', GM, GD, classifier, n, m, use_WKNKN, K, eta);

% Display evaluation results
save('Ranking matrix.mat', 'matrixs'); % Save as MAT file
csvwrite('Ranking matrix.csv', matrixs); % Save as CSV file
fprintf('\n Predictions saved to file Ranking matrix.mat and Ranking matrix.csv\n\n');
disp('--------------------------------------------------------------');



