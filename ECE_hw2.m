rxn1 = [ [1,  1, -1, -1,  0,  0,  0,  0,  0]; ...
         [0,  0,  0,  1, -1,  0,  0,  0,  0];...
         [0,  0,  0,  0,  0, -1,  1,  0,  0];...
         [0,  0,  1, -1,  0,  0,  0,  0,  0];...
         [0,  0,  0,  0,  1,  0, -1,  0, -1];...
         [0,  0,  0,  0,  0,  0,  0, -1,  1]; ];

rxn2 = [ [1,  1, -1, -1,  0,  0,  0,  0];...
         [0,  0,  0,  1, -1,  0,  0,  0];...
         [0,  0,  0,  0,  0, -1,  1,  0];...
         [0,  0,  1, -1,  0,  0,  0,  0];...
         [0,  0,  0,  0,  1,  0, -1,  0];...
         [0,  0,  0,  0,  0,  0,  0, -1]; ];
     
rk1 = rank(rxn1);
rk2 = rank(rxn2);
Z1 = null(rxn1,'r');
Z2 = null(rxn2,'r');
rz1 = rank(Z1);
rz2 = rank(Z2);

K = [ [1,  0];
      [0, -1];
      [1,  1];
      [1,  0]; ];
rK = rank(K);

% new basis for each system
B1 = Z1;
B1(:,1) = B1(:,2) + B1(:,1);
rb1 = rank(B1);
B2 = Z2;
B2(:,1) = B2(:,2) + B2(:,1);
rb2 = rank(B2);

% flux balance analysis with constraints
f1 = [0; 0; 0; 0; -1; 0; 1; -1; 0];
lb1 = [0.05; 0; 0; 0; 0; 0.1; 0; 0; 0];
A1 = [ [0, 1, 0, 0, 0, 0, 0, 0, 0];...
       [0, 0, 0, 1, 0, 0, 0, 0, 0]; ...
       [0, 0, 0, 0, 0, 1, 0, 0, 0]; ...
       [0, 0, 0, 0, 0, 0, 0, 1, 0]; ...
       [0, 0, 0, 0, 0, 0, 0, 0, 1]; ];
b1 = [1; 1; 0.2; 1; 1];
Aeq1 = [ [1, 0, 0, 0, 0, 0, 0, 0, 0];
         rxn1];
beq1 = [0.05; 0; 0; 0; 0; 0; 0];
[x1,sc1] = linprog(f1,A1,b1,Aeq1,beq1,lb1);


f2 = [0; 0; 0; 0; -1; 0; 1; -1];
lb2 = [0.05; 0; 0; 0; 0; 0.1; 0; 0];
A2 = [ [0, 1, 0, 0, 0, 0, 0, 0];...
       [0, 0, 0, 1, 0, 0, 0, 0]; ...
       [0, 0, 0, 0, 0, 1, 0, 0]; ...
       [0, 0, 0, 0, 0, 0, 0, 1]; ];
b2 = [1; 1; 0.2; 1];
Aeq2 = [ [1, 0, 0, 0, 0, 0, 0, 0];
         rxn2];
beq2 = [0.05; 0; 0; 0; 0; 0; 0];
[x2,sc2] = linprog(f2,A2,b2,Aeq2,beq2,lb2);