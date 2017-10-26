function element = shape_coef(mesh,element)
%This function compute the shape coeficient for each element

Ne=mesh.Ne;

% Node coordinates
x=mesh.node.x(element.num_node(:,:))*1000;
y=mesh.node.y(element.num_node(:,:))*1000;

% Jacobian
jacobian=x(:,2).*y(:,3)+x(:,3).*y(:,1)+x(:,1).*y(:,2)-x(:,2).*y(:,1)-x(:,3).*y(:,2)-x(:,1).*y(:,3);

% b1,2,3 et c1,2,3 are the shape coef.
b1=(y(:,2)-y(:,3))./jacobian; %b1
b2=(y(:,3)-y(:,1))./jacobian; %b2
b3=(y(:,1)-y(:,2))./jacobian; %b3
c1=(x(:,3)-x(:,2))./jacobian; %c1
c2=(x(:,1)-x(:,3))./jacobian; %c2
c3=(x(:,2)-x(:,1))./jacobian; %c3

% surface 
surface=abs(jacobian)./2;

% output
element.shape_coef=[b1,b2,b3,c1,c2,c3];
element.surf=surface;

element.x=mean(mesh.node.x(element.num_node(:,1:3)),2);
element.y=mean(mesh.node.y(element.num_node(:,1:3)),2);