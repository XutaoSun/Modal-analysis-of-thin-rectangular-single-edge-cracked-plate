%rectangular plate, all boundaries are simply supported, the central crack is parallel to the upper and lower boundaries. 2 segments are used.
clear all
clc
format long
symm=1; %symmetry indicator, 0: symmetric, 1: antisymmetric
n=38; %the number of terms of admissible functions in each of the two coordinate directions
n2=n^2;
a=2*0.32; 
b=0.32; %a, b are the plate dimensions in x and y directions
h=0.001*a; %thickness
c=0.4*a; %crack length
a1=a-c;
a2=c;
b1=b/2;
b2=b/2; %dimensions of each segment
al1=a1/b1;
al2=a2/b2; %aspect ratios of segments
nu=0.3; %the Poisson's ratio
E=69e9; %Young's modulus
Rho=2810; %density
D=E*h^3/(12*(1-nu^2)); %bending rigidity of the plate
%penalty terms for edges of the half plate, kp1 kr1(left edge) kp2 kr2(upper edge) kp3 kr3(right edge) kp4 kr4(lower edge)
kp1=1e9; %1e9; %left edge penalty
kp2=1e9; %upper edge
kp3=1e9; %1e9; %right edge
kr1=0; %left edge
kr2=0; %upper edge
kr3=0; %right edge
if symm==0 %symmetric
	ks4=1e9; 
	kp4=0;
	kr4=1e9; %lower edge, slip shear condition
else %antisymmetric
	ks4=0;
	kp4=1e9; 
	kr4=0; %lower edge, simply supported
end
%penalty terms for continuity of translation and rotation (kct and kcr) between segments
kct12=1e9;
kcr12=1e9;
%%%%%%%%%%%%%%%%%%%%%% define integrals %%%%%%%%%%%%%%%%%%%%
% Product of intergrals of admissible functions Emi00 stands for the product of two zeroth derivatives, Emi02 stands for the product of the zeroth derivative and second derivative (the order should be consistent with arguments in brackets)
Emi00(1,1)=1;
Emi00(1,2)=1/2;
Emi00(2,1)=1/2;
Emi00(3,1)=1/3;
Emi00(1,3)=1/3;
Emi00(2,2)=1/3;
Emi00(2,3)=1/4;
Emi00(3,2)=1/4;
Emi00(3,3)=1/5;
for i=4:n
    Emi00(2,i)=-(1-cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(i,2)=-(1-cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(3,i)=(2*cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(i,3)=(2*cos((i-3)*pi))/((i-3)^2*pi^2);
end
for i=4:n
    Emi00(i,i)=1/2;
end
%%%%
Emi22(3,3)=4;
for i=4:n
    Emi22(i,i)=0.5*(i-3)^4*pi^4;
end
%%%%
Emi02(1,3)=2;
Emi02(2,3)=1;
Emi02(3,3)=2/3;
for i=4:n
    Emi02(2,i)=1-cos((i-3)*pi);
    Emi02(3,i)=-2*cos((i-3)*pi);
    Emi02(i,i)=-0.5*(i-3)^2*pi^2;
end
%%%%
Emi11(2,2)=1;
Emi11(2,3)=1;
Emi11(3,2)=1;
Emi11(3,3)=4/3;
for i=4:n
    Emi11(2,i)=cos((i-3)*pi)-1;
    Emi11(i,2)=cos((i-3)*pi)-1;
    Emi11(3,i)=2*cos((i-3)*pi);
    Emi11(i,3)=2*cos((i-3)*pi);
    Emi11(i,i)=0.5*(i-3)^2*pi^2;
end
%%%%
Emi20=Emi02';
Fni00=Emi00;
Fni20=Emi20;
Fni02=Emi02;
Fni11=Emi11;
Fni22=Emi22;
%%%%
%%%%%%%%%%%%%%%%%% mass matrix %%%%%%%%%%%%%%%%%%
for m=1:n 
	for i=1:n 
		u=i+(m-1)*(n); 
		for nn=1:n
			for j=1:n
				v=j+(nn-1)*(n);
				M(u,v)=Emi00(m,nn)*Fni00(i,j)*a1*b1*Rho*h;
				M(u+n2,v+n2)=Emi00(m,nn)*Fni00(i,j)*a2*b2*Rho*h;
			end
		end
	end
end
%%%%%%%%%%%%%%%%%% define values of admissible functions at edges or crack %%%%%%%%%%%%%
%KP1是板的左边界，KP2是板的右边界，KP3是下边界，KP4是上边界，KC1是左边界转角，KC2是右边界转角，KC3是下边界转角，KC4是上边界转角。
%KS3是下边界slip shear condition
for i=1:n
	KP1(i)=1;
	KP3(i)=1;
	KC1(i)=0;
	KC2(i)=0;
	KS3(i)=0; %下边界slip shear condition
end
KP1(2)=0;
KP1(3)=0;
KP3(2)=0;
KP3(3)=0;
KC1(2)=1;
KC2(1)=0;
KC2(2)=1;
KC2(3)=2;
KC3=KC1;
KC4=KC2;
KP2(1)=1;
KP4(1)=1;
KP2(2)=1;
KP4(2)=1;
KP2(3)=1;
KP4(3)=1;
for i=4:n
	KP2(i)=cos((i-3)*pi);
	KP4(i)=KP2(i);
end
%%%%%%%%%%%%%%%%%%% stiffness matrix %%%%%%%%%%%%%%%%%%%
K=zeros(2*n^2,2*n^2);
for m=1:n
	for i=1:n
		u=i+(m-1)*n;
		for nn=1:n
			for j=1:n
				v=j+(nn-1)*n;
				%there are 2 plate segments, so imagine the dimension of the stiffness matrix K is 2X2
				%for K11
				K(u,v)=K(u,v)+(Emi22(m,nn)*Fni00(i,j)+al1^4*Emi00(m,nn)*Fni22(i,j)+nu*al1^2*Emi02(m,nn)*Fni20(i,j)+nu*al1^2*Emi20(m,nn)*Fni02(i,j)+2*(1-nu)*al1^2*Emi11(m,nn)*Fni11(i,j))*D*b1/a1^3; %the segment 1 itself
				K(u,v)=K(u,v)+kct12*KP2(m)*KP2(nn)*Fni00(i,j)*b1; %part of the translational boundary condition between segment 1 and segment 2 (right boundary of segment 1)
				K(u,v)=K(u,v)+kcr12*KC2(m)/a1*KC2(nn)/a1*Fni00(i,j)*b1; %part of the rotational boundary condition between segment 1 and segment 2
				
				K(u,v)=K(u,v)+kp1*KP1(m)*KP1(nn)*Emi00(i,j)*b1; %the translational boundary condition for the left edge of segment 1
				K(u,v)=K(u,v)+kr1*KC1(m)/a1*KC1(nn)/a1*Emi00(i,j)*b1; %the rotational boundary condition for the left edge of segment 1
				
				K(u,v)=K(u,v)+kp2*KP4(i)*KP4(j)*Emi00(m,nn)*a1; %the translational boundary condition for the upper edge of segment 1
				K(u,v)=K(u,v)+kr2*KC4(i)/b1*KC4(j)/b1*Emi00(m,nn)*a1; %the rotational boundary condition for the upper edge of segment 1
				K(u,v)=K(u,v)+kp4*KP3(i)*KP3(j)*Emi00(m,nn)*a1; %the translational boundary condition for the lower edge of segment 1
				K(u,v)=K(u,v)+ks4*KS3(i)*KS3(j)*Emi00(m,nn)*a1; %the slip shear boundary condition for the lower edge of segment 1
				K(u,v)=K(u,v)+kr4*KC3(i)/b1*KC3(j)/b1*Emi00(m,nn)*a1; %the rotational boundary condition for the lower edge of segment 1
				%for K22
				K(u+n2,v+n2)=K(u+n2,v+n2)+(Emi22(m,nn)*Fni00(i,j)+al2^4*Emi00(m,nn)*Fni22(i,j)+nu*al2^2*Emi02(m,nn)*Fni20(i,j)+nu*al2^2*Emi20(m,nn)*Fni02(i,j)+2*(1-nu)*al2^2*Emi11(m,nn)*Fni11(i,j))*D*b1/a2^3; %the segment 2 itself
				K(u+n2,v+n2)=K(u+n2,v+n2)+kct12*KP1(m)*KP1(nn)*Fni00(i,j)*b1; %part of the translational boundary condition between segment 2 and segment 1 (left boundary of segment 1)
				K(u+n2,v+n2)=K(u+n2,v+n2)+kcr12*KC1(m)/a2*KC1(nn)/a2*Fni00(i,j)*b1; %part of the rotational boundary condition between segment 2 and segment 1
				
				K(u+n2,v+n2)=K(u+n2,v+n2)+kp3*KP2(m)*KP2(nn)*Emi00(i,j)*b1; %the translational boundary condition for the right edge of segment 2
				K(u+n2,v+n2)=K(u+n2,v+n2)+kr3*KC2(m)/a2*KC2(nn)/a2*Emi00(i,j)*b1; %the rotational boundary condition for the right edge of segment 2
				
				K(u+n2,v+n2)=K(u+n2,v+n2)+kp2*KP4(i)*KP4(j)*Emi00(m,nn)*a2; %the translational boundary condition for the upper edge of segment 2
				K(u+n2,v+n2)=K(u+n2,v+n2)+kr2*KC4(i)/b1*KC4(j)/b1*Emi00(m,nn)*a2; %the rotational boundary condition for the upper edge of segment 2
				%for K12
				K(u,v+n2)=K(u,v+n2)-kct12*KP2(m)*KP1(nn)*Fni00(i,j)*b1; %part of the translational boundary condition between segment 1 and segment 2
				K(u,v+n2)=K(u,v+n2)-kcr12*KC2(m)/a1*KC1(nn)/a2*Fni00(i,j)*b1; %part of the rotational boundary condition between segment 1 and segment 2
				%for K21
				K(u+n2,v)=K(u+n2,v)-kct12*KP1(m)*KP2(nn)*Fni00(i,j)*b1; %part of the translational boundary condition between segment 1 and segment 2
				K(u+n2,v)=K(u+n2,v)-kcr12*KC2(nn)/a1*KC1(m)/a2*Fni00(i,j)*b1; %part of the rotational boundary condition between segment 1 and segment 2
			end
		end
	end
end
[X,eigenvalues] = eig(K,M);
eigenvalues = diag(eigenvalues);
Omega=eigenvalues.^0.5*a^2*(Rho*h/D)^0.5;
Omega=real(Omega)
x_nodes=linspace(0,a,31);
y_nodes=linspace(0,b/2,31); %coordinates of nodes in the half segment
X_nodes=linspace(0,a,31);
Y_nodes=linspace(-b/2,b/2,61); %coordinates of nodes in the full plate, keep the same coordinates as in the half segment but expand coordinates to the other half segment. It is crucial to use the same coordinate system and scale. 
for MODE_NUMBER=1:5	
	for i=1:size(x_nodes,2)
		x=x_nodes(i);
		if x<=a1 %when the point is in the left segment (segment 1)
			EIGENVECTOR=real(X(1:n2,MODE_NUMBER))';
			for j=1:size(y_nodes,2)
				y=y_nodes(j);
				SET_X(1)=1;
				SET_X(2)=x/a1;
				SET_X(3)=(x/a1)^2;
				SET_Y(1)=1;
				SET_Y(2)=y/b1;
				SET_Y(3)=(y/b1)^2;
				for k=1:n-3
					SET_X(k+3)=cos(k*pi*x/a1);
					SET_Y(k+3)=cos(k*pi*y/b1);
				end
				counter=0;
				for p=1:n
					for q=1:n
						counter=counter+1;
						EVALUATED_FUNCTIONS(counter)=SET_X(p)*SET_Y(q);
					end
				end
				EVALUATED_MODE(i,j)=-dot(EVALUATED_FUNCTIONS,EIGENVECTOR);
			end
		else %when the point is in the right segment (segment 2)
			EIGENVECTOR=real(X(n2+1:2*n2,MODE_NUMBER))';
			for j=1:size(y_nodes,2)
				y=y_nodes(j);
				SET_X(1)=1;
				SET_X(2)=(x-a1)/a2;
				SET_X(3)=((x-a1)/a2)^2;
				SET_Y(1)=1;
				SET_Y(2)=y/b1;
				SET_Y(3)=(y/b1)^2;
				for k=1:n-3
					SET_X(k+3)=cos(k*pi*(x-a1)/a2);
					SET_Y(k+3)=cos(k*pi*y/b1);
				end
				counter=0;
				for p=1:n
					for q=1:n
						counter=counter+1;
						EVALUATED_FUNCTIONS(counter)=SET_X(p)*SET_Y(q);
					end
				end
				EVALUATED_MODE(i,j)=-dot(EVALUATED_FUNCTIONS,EIGENVECTOR);
			end
		end
	end
	figure
	surf(x_nodes,y_nodes,EVALUATED_MODE) %mode shape of half plate
	if symm==0 %symmetric mode shape
		for i=1:size(X_nodes,2)
			for j=1:size(y_nodes,2)
				MODE(i,j)=EVALUATED_MODE(i,size(y_nodes,2)+1-j)';
			end
		end
		for i=1:size(X_nodes,2)
			for j=size(y_nodes,2)+1:size(Y_nodes,2)
				MODE(i,j)=EVALUATED_MODE(i,j-(size(y_nodes,2)-1))';
			end
		end
		figure
		Modeshape=surf(Y_nodes,X_nodes,MODE); %mode shape of full plate
		axis off;
		zdir=[0,90];
		rotate(Modeshape,zdir,180)
	else %antisymmetric mode shape
		for i=1:size(X_nodes,2)
			for j=1:size(y_nodes,2)
				MODE(i,j)=EVALUATED_MODE(i,size(y_nodes,2)+1-j)';
			end
		end
		for i=1:size(X_nodes,2)
			for j=size(y_nodes,2)+1:size(Y_nodes,2)
				MODE(i,j)=-EVALUATED_MODE(i,j-(size(y_nodes,2)-1))';
			end
		end
		figure
		Modeshape=surf(Y_nodes,X_nodes,MODE); %mode shape of full plate
		axis off;
		zdir=[0,90];
		rotate(Modeshape,zdir,180)
	end
end

%X_nodes=linspace(0,a,101);
%Y_nodes=linspace(-b/2,b/2,201);
%%symmetric mode shape
%for i=1:101
%	for j=1:101
%		MODE(i,j)=EVALUATED_MODE(i,j)';
%	end
%end
%for i=1:101
%	for j=102:201
%		MODE(i,j)=EVALUATED_MODE(i,j-100)';
%	end
%end
