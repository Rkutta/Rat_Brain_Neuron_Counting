% Script for Image Analysis
% Version 0.13
clear;
clc;
%% Load in Image and set parameters
input_image = 'LO 69 1-2 ACC Right.tif'; 
%input('Input image file name: ');
image = imread(input_image);
% threshold_upper = input('Input upper threshold for segmentation: ');
% threshold_lower = input('Input lower threshold for segmentation: ');

%% Convert to Grey Scale
imageG = rgb2gray(image);

%% Image Segmentation
mask = false(size(imageG)); % mask matrix for activecontour filter
[m,n] = size(imageG);
values = []; % matrix to store the pixel intensity values
i = 1;
% Get intensity values for pixels falling within the threshold
for x = 1:m
    for y = 1:n
        % Where the upper and lower thresholds come into play
        % for segmentation (
        if 0.01 < imageG(x,y) && imageG(x,y) < 70
            values(i) = imageG(x,y);
            i = i + 1;
        end
    end
end
% calculate the mean intensity of pixels withing the threshold range
mean_intensity = mean(values);
% fill mask funtion, flagging "true" pixel locations that are <= mean 
% intensity value

for x = 1:m
    for y = 1:n
        if imageG(x,y) < mean_intensity
            mask(x,y) = true;
        end
    end
end
%mask = imageG < mean_intensity;
% actual filtering using activecontour filter
BW = activecontour(imageG,mask);
% show filtered image
figure
imshow(BW)

%% Object Analysis
% Boundary analysis using bwboundaries function
% B = cell array of boundary pixel locations
% L = label matrix L where objects and holes are labeled
% N = number of objects found
% A - adjacency matrix
[B,L,N,A] = bwboundaries(BW,'noholes');
boundary_areas = []; % matrix for storing the areas of each boundary
for k = 1:length(B)
    boundary = B{k};
    % analyze area using polyarea function
    boundary_areas(k) = polyarea(boundary(:,2),boundary(:,1));
end
% from boundary_areas matrix, calculate average, and mean ares
average_area = mean(boundary_areas);
median_area = median(boundary_areas);
% calculated pixel area of a neuron
neuron_size = 598;
% upper and lower bounds for pixel area of nucleus of a neuron
lower_bound_nucleus_size = 0;
upper_bound_nucleus_size = 19.384;
% standared deviation of boundary areas
% S = std(boundary_areas);
% count variables to count valid boundaries after filtering by area
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
% arrays to store boundary pixel locations of valid boundaries
a_boundaries = []; % <= average area
m_boundaries = []; % <= median area
ns_boundaries = []; % <= pixel area of neuron
r_boundaries = []; % 0.538 pixels^2 <= boundary area <= 19.384 pixels^2
% filter boundaries by area's <= average area
for k = 1:length(B)
    boundary = B{k};
    if boundary_areas(k) <= (average_area)
        a_boundaries{k} = boundary;
        count1 = count1 + 1;
    end
end
% filter boundaries by area's <= median area
for k = 1:length(B)
    boundary = B{k};
    if boundary_areas(k) <= (median_area)
        m_boundaries{k} = boundary;
        count2 = count2 + 1;
    end
end
% filter boundaries by area's <= pixel area of neuron
for k = 1:length(B)
    boundary = B{k};
    if boundary_areas(k) <= (neuron_size)
        ns_boundaries{k} = boundary;
        count3 = count3 + 1;
    end
end
% filter boundaries by area's >= 0 pixels^2 and <= 19.384 pixels^2
for k = 1:length(B)
    boundary = B{k};
    if (lower_bound_nucleus_size <= boundary_areas(k)) && (boundary_areas(k) <= upper_bound_nucleus_size)
        r_boundaries{k} = boundary;
        count4 = count4 + 1;
    end
end
%% Object Plotting
% Plot objects outlines imposed on chosen image (currently set to greyscale
% image)
plots = [];
plots(1:5) = 0;
in = 'Blank';
% Display user instructions
fprintf('Object Plotting\n\n\')
disp('Enter 1 to plot objects filterd by average area')
disp('Enter 2 to plot objects filtered by median area')
disp('Enter 3 to plot objects filtered by pixel area of a neuron')
disp('Enter 4 to plot objects filtered by pixel area range for nucleus of a neuron')
disp('Enter 5 to plot all detected objects')
disp('Enter "all" to plot all object filterings (use single quotes)')
disp('Enter "plot" to make plots (use single quotes)')
disp('Enter "skip" to skip plotting (use single quotes)')
fprintf('\n')
% User input
while ~strcmp(in,'plot')
    in = input('Input command -> ');
    if in == 1
        plots(in) = 1;
    end
    if in == 2
        plots(in) = 1;
    end
    if in == 3
        plots(in) = 1;
    end
    if in == 4
        plots(in) = 1;
    end
    if in == 5
        plots(in) = 1;
    end
    if strcmp(in,'all')
        plots(1:5) = 1;
        break
    end
    if strcmp(in,'skip')
        break
    end
end  
% Plot objects filtered by average area
if plots(1)
    figure
    imshow(imageG)
    hold on
    title(['Dot Count Using Average Area: ',num2str(count1)])
    for k = 1:length(a_boundaries)
        boundary = a_boundaries{k};
        if ~isempty(boundary)
            plot(boundary(:,2),boundary(:,1),'Color','red','LineWidth',1)
        end
    end
end
% Plot objects filtered by median area
if plots(2)
    figure
    imshow(imageG)
    hold on
    title(['Dot Count Using Median Area: ',num2str(count2)])
    for k = 1:length(m_boundaries)
        boundary = m_boundaries{k};
        if ~isempty(boundary)
            plot(boundary(:,2),boundary(:,1),'Color','red','LineWidth',1)
        end
    end
end
% Plot objects filtered by pixel area of a neuron
if plots(3)
    figure
    imshow(imageG)
    hold on
    title(['Dot Count Using Calculated Pixel Area of Neruon: ',num2str(count3)])
    for k = 1:length(ns_boundaries)
        boundary = ns_boundaries{k};
        if ~isempty(boundary)
            plot(boundary(:,2),boundary(:,1),'Color','red','LineWidth',1)
        end
    end
end
% Plot objects filtered by nucleus area of a neuron
if plots(4)
    figure
    imshow(imageG)
    hold on
    title(['Dot Count Using Calculated Pixel Range for Nucleus of Neuron: ', num2str(count4)])
    for k = 1:length(r_boundaries)
        boundary = r_boundaries{k};
        if ~isempty(boundary)
            plot(boundary(:,2),boundary(:,1),'Color','red','LineWidth',1)
        end
    end
end
% Plot all objects
if plots(5)
    figure
    imshow(imageG)
    hold on
    title(['Dot Count Using All Areas: ',num2str(N)])
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'Color','red','LineWidth',1)
    end
end
%% Data Storage for Number of Objects
% count1 = filtered by average object area
% count2 = filtered by mean object area
% count3 = filtered by neuron size
% count4 = filtered by neuron nucleus size
% N =  unfiltered
counts = [count1 count2 count3 count4 N];