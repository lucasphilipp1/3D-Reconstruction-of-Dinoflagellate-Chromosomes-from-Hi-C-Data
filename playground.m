% get_model = importdata('test.xyz');
% get_model(:,1,:) = [];
get_model = cell2mat(struct2cell(load('matrix.mat')));

numPoints = size(get_model,1);
MyColor = linspace(1,numPoints,numPoints)';
%create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';

%plot the figure 
f=figure
hold on
plot3(get_model(:,1),get_model(:,2),get_model(:,3),'Color', [.6 .6 .6])
colormap jet
axis equal


save_H_values = zeros(1,2);
%number of monomers that specify each disc or spiral
skip = 36;
%ensures that points A,B,C,D are obtained before calculating H
count = 1;
%tracks the midpoints of AB and CD in order to calculate EF
track_midpoints = {};
%saves the vectors AB and CD to calculate H 
track_vectors = {};

%iterates through the structure to plot vectors and calculate H values 
for i=1:skip+(1/2*skip):numPoints
    i
    %gets points on the structure  
    point1 = get_model(i,:,:);
    point2 = get_model(i+(skip-1),:,:);
    %saves the midpoint between the two points
    track_midpoints{end+1} = (point1+point2)/2;
    %saves the vector that connects the two points
    track_vectors{end+1} = point2-point1;
    %plots the vector onto the structure for visualization 
    plot3([point1(1) point2(1)],[point1(2) point2(2)],[point1(3) point2(3)],'r-^', 'LineWidth',0.5);
    count = count + 1;
    %if two vectors are obtained, get EF and calculate H 
    if size(track_vectors, 2) == 2
        %get the vector EF given vectors AB and CD
        E = track_midpoints{1};
        F = track_midpoints{2};
        vector_EF = F-E;
        plot3([E(1) F(1)],[E(2) F(2)],[E(3) F(3)],'r-^', 'LineWidth',0.5);
        %calculate H
        calculate_H = (dot(vector_EF, cross(track_vectors{2}, track_vectors{1})))/(norm(track_vectors{1})*norm(track_vectors{2})*norm(vector_EF));
        save_H_values = cat(1,save_H_values,[i calculate_H]);
        %reset variables
        count=1; 
        track_midpoints = track_midpoints(2);
        track_vectors = track_vectors(2);
    end
end

save_H_values(1,:) = [];
%plots H value and the number of monomers iterated
hold off
plot(save_H_values(:,1),save_H_values(:,2))