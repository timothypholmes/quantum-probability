%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Timothy Holmes           
%                         quantum-probability
%                               1/1/18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Examples and inputs

%Uncomment for custom input
operator = [ 0 -1i 0 ; 1i 0 -1i; 0 1i 0 ]/sqrt(2);
state = [ 3 ; 1i ; -1i ];

%operators are Sz
%eigenvectors |+>
%eigenvalue hbar/2
%Sz = [1 0; 0 -1; 0 0];
%Sx = [0 1; 1 0];
%Sy = [0 -1i, 1i 0];
% plusEigVal = [1; 0];
% plusEigValX = (1/sqrt(2))*[1; 1];
% plusEigValY = (1/sqrt(2))*[1; 1i];
% minusEigVal = [0; 1];
% minusEigValX = (1/sqrt(2))*[1; -1];
% minusEigValY = (1/sqrt(2))*[1; -1i];

function quantum_probability_func(operator, state)
%% Defined inputes
A = operator; % Redefined inputes for simplicity 
B = state;

%% Normailze 
Norm = sqrt(1/dot(B,B));    %Normalize the state, find C
NormState = Norm.*state;    %Normalize the state

%% Results Calculation
e = eig(A); %Find the eigenvectors
Results = e; %Results of the operator

%% States Calculation
[V,D] = eig(A);
States = V; %The states of the eigenvector
Diagonalize = D; %The diangonalize 

%% Expected Value Calculation
Expected = (NormState' * A) * NormState; % calculation of <A^2>
ExpectedSq = (NormState' * (A^2)) * NormState; %calculation of <A>^2

%% Uncertainty Calculation
DeltaA = sqrt(ExpectedSq - (Expected^2)); %calculation sqrt( <A^2> - <A>^2)

%% Probability Calculation
Probs = abs((V'*NormState).^2); %Find the probability
    
%% Print Values (outputs of all calculations)
fprintf('The probabilities are; \n');
fprintf([repmat('%f\t', 1, size(Probs, 2)) '\n'],Probs);
fprintf('The resultes are; \n');
fprintf([repmat('%f\t', 1, size(Results, 2)) '\n'],Results);
fprintf('The states are(Top is real,bottom is compex); \n');
fprintf([repmat('%f\t', 1, size(States, 2)) '\n'],[real(States.'),imag(States.')]);
%fprintf([repmat('%f\t', 1, size(imag(States), 2)) '\n'],imag(States.'));
fprintf('The expected value is %1.3f. \n',Expected);
fprintf('The uncertainty is %1.3f. \n',DeltaA);
fprintf('The Diagonalized matrix is \n');
fprintf([repmat('%f\t', 1, size(Diagonalize, 2)) '\n'],Diagonalize);

end