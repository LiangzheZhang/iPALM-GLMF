 function y3=alg_iPALMGLMF(Y,GM,GD,use_WKNKN,K,eta)
%
% INPUT:
%  Y:           interaction matrix
%  GM:          microbe similarity matrix  
%  GD:          disease similarity matrix  
%
% OUTPUT:
%  y3:  prediction matrix
 
    %best parameters
    k = 50;
    lambda_l = 2^-2;
    lambda_d = 10^-1;
    lambda_t  = 10^-1;
    lambda_s  = 10^-2;
    
    %nearest neighbor graph (for sparsification of similarity matrices)
    p = 5;

    % stopping criterion (number of iterations)
    num_iter = 2;

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,GM,GD,K,eta);
    end

    % preprocessing Sd & St
    GM = preprocess_PNN(GM,p);
    GD = preprocess_PNN(GD,p);

    % Laplacian Matrices
    Dd = diag(sum(GM));
    Dt = diag(sum(GD));
    Ld = Dd - GM;
    Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));
    Lt = Dt - GD;
    Lt = (Dt^(-0.5))*Lt*(Dt^(-0.5));

    % initialize A & B
    [A,B] = initializer(Y,k);

    % obtain the feature matrix
    [A,B] = alg_iPALMGLMF_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,lambda_s,num_iter);

    % compute prediction matrix
    y3 = A*B';
    
    %--------------------------------------------------------------------

end