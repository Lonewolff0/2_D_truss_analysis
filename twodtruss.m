clc;
close all;
Q=0;
W=0;
%importing node data from excel
while Q==0
    fprintf('\nKEY-1-->TO VIEW FORMAT OF INPUT FILE\nKEY-2-->TO SELECT INPUT FILE\nANY OTHER KEY TO CANCEL\n')
    fprintf('____________________________________\n')
    prompt='INPUT KEY-';
    W = input(prompt);
    if W==1
        disp('_________________________________')
        disp('INPUT FILE SHOULD BE IN FORMAT OF xlsx WITH TWO SHEET-1,SHEET-2')
        disp('___________________________')
        disp('SHEET-1 SHOULD CONTAIN DATA ABOUT NODES')
        disp('FORMAT OF SHEET-1-')
        disp('______________________')
        disp('NM  X  Y  Fx  Fy  Mx  My')
        fprintf('NM-->Node Number\nX-->X COORDINATE OF NODE\nY-->Y COORDINATE OF NODE\nX-->X COMPONENT OF FORCE ACTING ON NODE\nY-->Y COMPONENT OF FORCE ACTING ON NODE\nMx = 1 IF NODE CAN MOVE IN X DIRECTION OR Mx=0\nMy = 1 IF NODE CAN MOVE IN Y DIRECTION OR My=0\n')
        disp('_________________________________')
        disp('SHEET-2 SHOULD CONTAIN DATA ABOUT ELEMENTS')
        disp('____________________')
        disp('FORMAT OF SHEET-2-')
        disp('______________________')
        disp('NODE_i  NODE_j  E  A')
        fprintf('Node_i-->ONE OF BOUNDARY NODE OF ELEMENT\nNODE_j-->OTHER BOUNDARY NODE OF ELEMENT\nE-->YOUNGS MODULUS OF ELEMENT\nA-->AREA OF CROSS SECTION OF ELEMENT\n')
        disp('_________________________________')
    elseif W==2
        [file,path] = uigetfile('*.xlsx');
        if isequal(file,0)
           disp('USER SELECTED CANCEL');   
           disp('--------------------------------------------------------')     
           disp('DO YOU WISH TO SELECT AGAIN THEN INPUT-1')
           disp('OR INPUT ANYKEY TO CANCEL')
           prompt='INPUT KEY-';
           inp = input(prompt);
           if inp==1
               Q=0;     
           else
               Q=10;    
           end
        else
           disp(['USER SELECTED-', fullfile(path,file)]);

            node_data = xlsread(file,1);
            element_data= xlsread(file,2);
            o=0; 
            node_number=sort(node_data(:,1));
            number_of_nodes=length(node_number);
            size_of_element_matrix = size(element_data);
            number_of_elements = size_of_element_matrix(1);
            element_boundary_nodes=element_data(:,1:2);
            for i=1:1:number_of_nodes
                if sum(node_number==node_number(i))~=1
                    o=1;
                    break
                end
                if o==1
                    break
                end
            end
            for i=1:number_of_elements-1
                boundary_nodes_of_element = element_boundary_nodes(i,:);
                for j=i+1:number_of_elements
                    boundary_nodes_of_element1 = element_boundary_nodes(j,:);
                    if (boundary_nodes_of_element1(1)==boundary_nodes_of_element(1) && boundary_nodes_of_element(2)==boundary_nodes_of_element1(2)) || (boundary_nodes_of_element1(1)==boundary_nodes_of_element(2) && boundary_nodes_of_element(1)==boundary_nodes_of_element1(2))
                        o=3;
                        break
                    end
                end
            end
            nc=node_data(:,2:3);
            for i=1:1:number_of_nodes-1
                ncd=nc(i,:);
                for j=i+1:1:number_of_nodes
                    ncd1=nc(j,:);
                    if ncd(1)==ncd1(1) && ncd(2)==ncd1(2)
                        o=2;
                        break
                    end
                end
                if o==2
                    break 
                end
            end
            if o==0
                node_number1=node_data(:,1);
                node_number=sort(node_data(:,1));
                node_data1=node_data;
                for i=1:1:number_of_nodes
                    h=node_data(:,1);
                    d=find(h==node_number(i));
                    node_data1(i,1:7)=node_data(d,:);
                end
                node_data=node_data1;
                x = node_data(:,2); y = node_data(:,3); Fx = node_data(:,4); Fy = node_data(:,5); Mx = node_data(:,6); My = node_data(:,7); 
                node_i = element_data(:,1); node_j = element_data(:,2); E = element_data(:,3); A = element_data(:,4);
                %finding length, angles and global stiffness matrix
                Global_stiffness_matrix = zeros(number_of_nodes*2);
                G_s_m1 =Global_stiffness_matrix;
                EAL =zeros(number_of_elements,1);
                l =zeros(number_of_elements,1);
                theta =zeros(number_of_elements,1);
                for i =1:number_of_elements
                    d=find(node_number==node_i(i));
                    f=find(node_number==node_j(i));
                    theta(i) = atand((y(f)-y(d))/(x(f)-x(d))); 
                    l(i) = sqrt((x(d)-x(f))^2 + (y(f)-y(d))^2);
                    EAL(i) = E(i)*A(i)/l(i);
                    c = cosd(theta(i));
                    s = sind(theta(i));
                    c2 = c^2;
                    s2 = s^2;
                    k =EAL(i)* [c2 c*s -c2 -c*s;c*s s2 -c*s -s2;-c2 -c*s c2 c*s;-c*s -s2 c*s s2];
                    p=[2*node_i(i)-1 2*node_i(i) 2*node_j(i)-1 2*node_j(i)];
                    K{i} = k;
                    G_s_m1(p,p)=k;
                    Global_stiffness_matrix = Global_stiffness_matrix+G_s_m1;
                    G_s_m1=zeros(number_of_nodes*2);
                end
                %FORCE MATRIX
                Force_matrix=zeros(8,1);
                for i=1:number_of_nodes
                    Force_matrix(2*i-1,1) = Fx(i);
                    Force_matrix(2*i,1) = Fy(i);
                end
                u=[];
                t=1;
                for i=1:number_of_nodes
                    if Mx(i)==1
                        u(t)=2*i-1;
                        t=t+1;
                    end
                    if My(i)==1
                        u(t)=2*i;
                        t=t+1;
                    end
                end
                G = Global_stiffness_matrix(u,u)\Force_matrix(u);
                U=zeros(number_of_nodes*2,1);
                for i =1:length(u)
                    U(u(i)) = G(i,1);
                end
                %Calculating reaction forces
                Reaction_Forces_matrix = Global_stiffness_matrix*U;
                Force_in_each_element_matrix = zeros(number_of_elements,1);
                for i=1:number_of_elements
                    p=[2*node_i(i)-1 2*node_i(i) 2*node_j(i)-1 2*node_j(i)];
                    EF = K{i}*U(p);
                    Force_in_each_element_matrix(i)= (EF(1)^2 + EF(4)^2)^(1/2);
                    element_data(i,5)=Force_in_each_element_matrix(i);
                    element_data(i,6)=Force_in_each_element_matrix(i)/A(i);
                end
                %new coordinates
                new_x=zeros(number_of_nodes);
                new_y=zeros(number_of_nodes);
                for i=1:number_of_nodes
                    new_x(i) = x(i)+(U(2*i-1));
                    new_y(i) = y(i)+(U(2*i));
                end
                for i=1:number_of_nodes
                    node_data(i,8)=U(2*i-1);
                    node_data(i,9)=U(2*i);
                    node_data(i,10)=Reaction_Forces_matrix(2*i-1);
                    node_data(i,11)=Reaction_Forces_matrix(2*i);
                end
            
                Z = 0;
                while Z==0
                    disp('INPUT KEY-1--> TO VIEW GRAPICAL REPRESENTATION OF TRUSS DEFORMATION')
                    disp('      KEY-2--> TO VIEW DEFORMATION OF EACH NODE (NODE_NUMBER  Ux  Uy)')
                    disp('      KEY-3--> TO VIEW FORCE MATRIX WITH RESPECT TO NODES(NODE_NUMBER  Fx  Fy)')
                    disp('      KEY-4--> TO VIEW FORCES IN EACH ELEMENT(NODEi NODEj FORCE)')
                    disp('      KEY-5--> TO VIEW STRESS IN EACH ELEMENT(NODEi NODEj STRESS)')
                    disp('      ANY OTHER KEY CANCEL')
                    disp('_________________________________')
                    prompt='INPUT KEY-';
                    disp('_________________________________')
                    Z = input(prompt);
                    if Z ==1
                        MaxArea=min(A);
                        for i=1:number_of_elements
                            hold on;
                            fig=plot(x([node_i(i),node_j(i)]),y([node_i(i),node_j(i)]),'b*-',new_x([node_i(i),node_j(i)]),new_y([node_i(i),node_j(i)]),'r*-');
                            d=find(node_number==node_i(i));
                            f=find(node_number==node_j(i));
                            text(x(node_i(i)),y(node_i(i)), num2str(node_number(d)),'Fontsize', 15)
                            text(x(node_j(i)),y(node_j(i)), num2str(node_number(f)),'Fontsize', 15)
                            text(new_x(node_i(i)),new_y(node_i(i)), num2str(node_number(d)),'Fontsize', 15)
                            text(new_x(node_j(i)),new_y(node_j(i)), num2str(node_number(f)),'Fontsize', 15)
                            thickness=2*A(i)/MaxArea;
                            set(fig,{'LineWidth'},{thickness;thickness})
                            set(fig,{'MarkerSize'},{7})
                            legend('BEFORE DEFORMATION','AFTER DEFORMATION')
                        end
                        Z=0;
                    elseif Z==2
                        disp('DEFORMATION IN X AND Y DIRECTION OF EAXH NODE')
                        disp('  NODE_NUMBER    Ux    Uy')
                        disp(node_data(:,[1,8,9]))
                        disp('_________________________________')
                        Z=0;
                    elseif Z==3
                        disp('FORCE ON EACH NODE IN X AND Y DIRECTION')
                        disp('  NODE_NUMBER    Ux     Uy')
                        disp(node_data(:,[1,10,11]))
                        disp('_________________________________')
                        Z=0;
                    elseif Z==4
                        disp('FORCE IN EACH ELEMENT')
                        disp('  Node_i    Node_j     Force')
                        disp(element_data(:,[1,2,5]))
                        disp('_________________________________')
                        Z=0;
                    elseif Z==5
                        disp('STRESS IN EACH ELEMENT')
                        disp('Node_i    Node_j   STRESS')
                        disp(element_data(:,[1,2,6]))
                        disp('_________________________________')
                        Z=0;
                    else
                        Z=-1;
                    end
                end
                disp('DO YOU WISH TO GO TO FIRST MENU THEN INPUT-1')
                disp('OR INPUT ANYKEY TO CANCEL')
                prompt='INPUT KEY-';
                inp = input(prompt);

                if inp==1
                    Q=0;
                else
                    Q=10;
                    disp('USER SELECTED CANCEL')
                end
            elseif o==1
                disp('--------------------------------------------------------------------------')
                disp('CANNOT COMPUTE....TWO OR MORE NODES HAVE SAME NUMBER!!!')
                disp('--------------------------------------------------------------------------')
                disp('DO YOU WISH TO GO TO FIRST MENU THEN INPUT-1')
                disp('OR INPUT ANY OTHER KEY TO CANCEL')
                prompt='INPUT KEY-';
                inp = input(prompt);

                if inp==1
                    Q=0;
                else
                    Q=10;
                    disp('USER SELECTED CANCEL')
                end
            elseif o==2
                disp('--------------------------------------------------------------------------')
                disp('CANNOT COMPUTE....TWO OR MORE NODES HAVE SAME COORDINATES!!!')
                disp('--------------------------------------------------------------------------')
                disp('DO YOU WISH TO GO TO FIRST MENU THEN INPUT-1')
                ddisp('OR INPUT ANY OTHER KEY TO CANCEL')
                prompt='INPUT KEY-';
                inp = input(prompt);

                if inp==1
                    Q=0;
                else
                    disp('USER SELECTED CANCEL')
                    Q=10;
                end
            elseif o==3
                disp('--------------------------------------------------------------------------')
                disp('CANNOT COMPUTE....TWO OR MORE ELEMENTS HAVE HAVE BOUNDARY NODES!!!')
                disp('--------------------------------------------------------------------------')
                disp('DO YOU WISH TO GO TO FIRST MENU THEN INPUT-1')
                disp('OR INPUT ANY OTHER KEY TO CANCEL')
                prompt='INPUT KEY-';
                inp = input(prompt);

                if inp==1
                    Q=0;
                else
                    disp('USER SELECTED CANCEL')
                    Q=10;
                end
            end
        end
    else
        disp('USER SELECTED CANCEL')
        break
    end
end
disp('-----------------------END-------------------------')