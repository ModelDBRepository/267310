%AUTHOR: Lisa Blum Moyse
%        lisa.blum-moyse@inria.fr
%
% REFERENCE: Blum Moyse & Berry. Modelling the modulation of cortical Up-Down state switching by astrocytes
%
% LICENSE: CC0 1.0 Universal

clear all; close all;
ini = -0.01;

for sig = [1.5]
    folder = '/home/lisa/Bureau/ubuntu2/stage_these/Article_M2/scriptsUPDOWN/spiking_model/data/idx_FP_Jai_sup_0_2000A/';
    %folder = '/home/lisa/Bureau/ubuntu2/stage_these/Article_M2/scriptsUPDOWN/spiking_model/data/idx_FP_Jai0_2000A/';
    
    %adapt the grid
    [xq,yq] = meshgrid(ini:1*10^(-5):0.015, ini:1*10^(-5):0.015);

    T = readtable(char(strcat(folder,'Ce',num2str(round(sig*1000,1)),'.csv')));

    xe = table2array(T(:,1));
    ve = table2array(T(:,2));
    ye = table2array(T(:,3));
    
    T = readtable(char(strcat(folder,'Ci',num2str(round(sig*1000,1)),'.csv')));
    xi = table2array(T(:,1));
    vi = table2array(T(:,2));
    yi = table2array(T(:,3));

    T = readtable(char(strcat(folder,'Ca',num2str(round(sig*1000,1)),'.csv')));
    va = table2array(T(:,1));
    xa = table2array(T(:,2));
    ya = table2array(T(:,3));

%     figure();
%     scatter3(va,xa,ya,20,'.g');
%     hold on;
%     scatter3(xe,ve,ye,20,'.r');
%     scatter3(xi,vi,yi,20,'.b');
%     xlabel('re');
%     ylabel('ri');
%     zlabel('ra');
%     title('sig=',sig);

    vqe = griddata(xe,ye,ve,xq,yq,'natural');
    vqi = griddata(xi,yi,vi,xq,yq,'natural');
    vqa = griddata(xa,ya,va,xq,yq,'natural');

    figure();
    mesh(xq,vqe,yq,'edgecolor','r');
    hold on;
    mesh(xq,vqi,yq,'edgecolor','b');
    %scatter3(va,xa,ya,10,'g');
    mesh(xq,vqa,yq,'edgecolor','g');
    title('sig=',sig);
    xlabel('re');
    ylabel('ri');
    zlabel('ra');

end



