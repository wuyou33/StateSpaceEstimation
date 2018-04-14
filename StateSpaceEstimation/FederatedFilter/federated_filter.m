function [ globalState, globalCov, state1New, cov1New, state2New, cov2New ] = federated_filter( state1, cov1, state2, cov2 )
    % federated_filter. Federated filter.
    %   Utilised for fusing the data from two separate filters using folliwng expression
    %
    %   [ globalState, globalCov, state1New, cov1New, state2New, cov2New ] = federated_filter( state1, cov1, state2, cov2 )
    %
    %   Xg = [(P_1)^-1 + (P_2)^-1]^-1 * ( (P_1)^-1*X_1 + (P_2)^-1*X_2 );
    %   Pg = [(P_1)^-1 + (P_2)^-1]^-1;
    %
    %   where
    %       Xg   global state;
    %       Pg   global covariance;
    %       X_1  state estimation from 1-st local filter;
    %       P_1  covariance estimation from 1-st local filter;
    %       X_2  state estimation from 2-nd local filter;
    %       P_2  covariance estimation from 2-nd local filter.
    %
    %   The global estimation is fed back to two local filters by the following equations:
    %   Xi = Xg;
    %   Pi = (bi * (Pg)^-1)^-1;
    %
    %               ||Pi||^-1
    %   bi = ------------------------
    %             Sum(||Pi||^-1)
    %
    %   where
    %       i - number of local filter;
    %       ||*|| - Frobenius-norm.
    %
    %   INPUT
    %       state1 - state estimation from 1-st local filter;
    %       cov1   - covariance estimation from 1-st local filter;
    %       state2 - state estimation from 2-nd local filter;
    %       cov2   - covariance estimation from 2-nd local filter.
    %
    %   OUTPUT
    %       globalState - global state;
    %       globalCov   - global covariance;
    %       state1New   - corrected state for 1-st local filter;
    %       cov1New     - corrected covariance for 1-st local filter;
    %       state2New   - corrected state for 2-nd local filter;
    %       cov2New     - corrected covariance for 2-nd local filter.
    %
    narginchk(4, 4);
    
    if size(state1) ~= size(state2)
        error('[ federated_filter ] dimension mismatch between state1 and state2');
    end
    
    if size(cov1) ~= size(cov2)
        error('[ federated_filter ] dimension mismatch between cov1 and cov2');
    end
    
    invCov1 = inv(cov1);
    invCov2 = inv(cov2);
    
    globalState = (invCov1 + invCov2) \ (cov1 \ state1 + cov2 \ state2);
    globalCov   = inv(invCov1 + invCov2);
    
    state1New = globalState;
    state2New = globalState;
    
    invGlobalCov = inv(globalCov);
    invFrobeniusNorm1 = inv(norm(cov1, 'fro'));
    invFrobeniusNorm2 = inv(norm(cov2, 'fro'));
    
    cov1New = (invFrobeniusNorm1 \ (invFrobeniusNorm1 + invFrobeniusNorm2)) \ invGlobalCov;
    cov2New = (invFrobeniusNorm2 \ (invFrobeniusNorm1 + invFrobeniusNorm2)) \ invGlobalCov;
end
