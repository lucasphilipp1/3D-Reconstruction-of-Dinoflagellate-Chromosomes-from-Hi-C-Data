%get_model = importdata('test.xyz');
%get_model(:,1,:) = [];
get_model = cell2mat(struct2cell(load('matrix.mat')));

numPoints = size(get_model,1);
MyColor = linspace(1,numPoints,numPoints)';
% create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';
save_H_values = zeros(1,2);

skip=37; %each spiral 

f=figure
hold on
plot3(get_model(:,1),get_model(:,2),get_model(:,3),'Color', [.6 .6 .6])
colormap jet
axis equal
% xlim([min(get_model(:,1))*2 max(get_model(:,1))*2])
% ylim([min(get_model(:,2))*2 max(get_model(:,2))*2])
% zlim([min(get_model(:,3))*2 max(get_model(:,3))*2])
skip = 90;
count = 1;
track_midpoints = {};
track_vectors = {};
num_monomers = 0;
for j=1:1:numPoints-((skip+(1/2*skip))*2)
    for i=j:(skip+(1/2*skip)):j+((skip+(1/2*skip))*2)-1
        point1 = get_model(i,:,:);
        point2 = get_model(i+(skip-1),:,:);
        track_midpoints{end+1} = (point1+point2)/2;
        track_vectors{end+1} = point2-point1;
        plot3([point1(1) point2(1)],[point1(2) point2(2)],[point1(3) point2(3)],'r-^', 'LineWidth',0.5);
        count = count + 1;
        if count == 3
            %calculate H
            E = track_midpoints{1};
            F = track_midpoints{2};
            vector_EF = F-E;
            plot3([E(1) F(1)],[E(2) F(2)],[E(3) F(3)],'r-^', 'LineWidth',0.5);
            calculate_H = (dot(vector_EF, cross(track_vectors{2}, track_vectors{1})))/(norm(track_vectors{1})*norm(track_vectors{2})*norm(vector_EF));
            %reset variables
            count=1; 
            track_midpoints = {};
            track_vectors = {};
        end
    end
    save_H_values = cat(1,save_H_values,[j calculate_H]);
end
set = 0;
save_H_values(1,:) = [];
% hold off
% plot(save_H_values(:,1),save_H_values(:,2))
 