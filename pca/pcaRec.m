%%this is a program using pca for dimension reduction
%%x_data contains the data after dimension reduction
%%y is the matrix to reduce the dimension

[y,z,yy]=princomp(x);
y=y(:,1:12);
x=x';
[m,n]=size(x);
x_mean=mean(x,2);
x_var=(x-repmat(x_mean,1,n));
x_data=y'*x_var;
x_data=x_data';