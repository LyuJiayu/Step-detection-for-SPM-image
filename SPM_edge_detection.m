clc
clear
% Input
files = "Ar30\" + ["1.asc","2.asc","3.asc","4.asc","5.asc"];
% Image size
xysize = zeros(5,1);
%%

ans_1vertical = zeros(5,1);
ans_2horizontal = zeros(5,1);
ans_3_v_h =  zeros(5,1);
ans_4sum = zeros(5,1);
ans_11vertical = zeros(5,1);
ans_22horizontal = zeros(5,1);
ans_33_v_h =  zeros(5,1);
ans_44sum = zeros(5,1);
ans_5lossR = zeros(5,1);
ratio_up = zeros(5,1);
ratio_down = zeros(5,1);

threshold = {[0.05 0.12],[0.05 0.12],[0.05 0.12],[0.05 0.12],[0.05 0.12]};

prT_name = cell(5,1);

for h = 1:length(files)
    filename = files(h);
    prT_name(h) = get_file_name(files(h));
    I_raw = readmatrix(filename,"FileType","text",'NumHeaderLines',18);
    I_raw = imresize(I_raw,[512,512]); 
    I_raw = imgaussfilt(I_raw,2);
    I_rescale = rescale(I_raw); 
    xysize(h) = get_pic_length(filename); 
    [...
    ans_1vertical(h),...
    ans_2horizontal(h),...
    ans_3_v_h(h),...
    ans_4sum(h),...
    ans_11vertical(h),...
    ans_22horizontal(h),...
    ans_33_v_h(h),...
    ans_44sum(h),...
    ans_5lossR(h),...
    L,...
    L_1vertical,...
    L_2horizontal,...
    P...
    ] = process_SPM(h,I_raw,I_rescale,xysize(h),threshold{h});
    
end

% Plot Result
tiledlayout(1,5);
for i = 1:5
    nexttile
    p = imread("f"+num2str(i)+".png");
    imshow(p)
    title({...
        prT_name(i)+"    "+ num2str(xysize(i)) + " nm"...
        })
    
    xlabel({...
        "Threshold = [" + num2str(threshold{i}(1)) + "," + num2str(threshold{i}(2)) + "]";...
        sprintf("Step Density = %.4f nm/nm^-^2",ans_11vertical(i)+ans_22horizontal(i));...
        sprintf("V_a_t_o_m = %.4f 10^1^5cm^-^2",ans_11vertical(i)/3.25);...
        sprintf("H_a_t_o_m = %.4f 10^1^5cm^-^2",ans_22horizontal(i)/5.25);...
        sprintf("V_a_t_o_m/H_a_t_o_m = %.2f",ans_11vertical(i)/3.25/(ans_22horizontal(i)/5.25));...
        },...
        'Position',[400 820]...
        )

end
% Save as excel
output_excel = zeros(4,1);
for i = 1:length(ans_1vertical)
    output_excel(4*i-3) = ans_1vertical(i);
    output_excel(4*i-3+1) = ans_2horizontal(i);
    output_excel(4*i-3+2) = ans_3_v_h(i);
    output_excel(4*i-3+3) = ans_4sum(i);
end
output_excel_2 = zeros(4,2);
for i = 1:length(ans_11vertical)
    output_excel_2(4*i-3,2) = ans_11vertical(i);
    output_excel_2(4*i-3+1,2) = ans_22horizontal(i);
    output_excel_2(4*i-3+2,2) = ans_33_v_h(i);
    output_excel_2(4*i-3+3,2) = ans_44sum(i);
    
    output_excel_2(4*i-3,1) = xysize(i)^2;
end
%%
% Output
% Calculate the weighted average
[mean_StepDensity,error_S] = weighted_mean(ans_4sum,xysize.^2)
% vertical step atom number
[mean_V_atoms,error_V_a] = weighted_mean(ans_11vertical*0.1/0.325,xysize.^2) % 0.325 nm for <1-210> step
% horizontal step atom number
[mean_H_atoms,error_H_a] = weighted_mean(ans_22horizontal*0.1/0.525,xysize.^2)  % 0.525 nm for <0001> step
% total atom number
[mean_T_a,error_T_a] = weighted_mean(ans_11vertical*0.1/0.325+ans_22horizontal*0.1/0.525,xysize.^2)
%%
function [...
    ans_1vertical,...
    ans_2horizontal,...
    ans_3_v_h,...
    ans_4sum,...
    ans_11vertical,...
    ans_22horizontal,...
    ans_33_v_h,...
    ans_44sum,...
    ans_5lossR,...
    L,...
    L_1vertical,...
    L_2horizontal,...
    P] = process_SPM(h,I_raw,I_rescale,xysize,threshold)

[E,~] = edge(I_rescale,"Canny",threshold); % Edge detection

result1 = zeros([512,512]);
result2 = zeros([512,512]);
for row = 1:512
    for col = 1:512
        if E(row,col) == 1
            %             row,col;
            [a1,pos1] = findbound(E,row,col,1);
            [a2,pos2] = findbound(E,row,col,2);
            [former1,latter1] = findnei(a1,pos1);
            [former2,latter2] = findnei(a2,pos2);
            n1 = round((max(I_raw(row-former1:row+latter1,col))-min(I_raw(row-former1:row+latter1,col)))/0.3);
            n2 = round((max(I_raw(row,col-former2:col+latter2))-min(I_raw(row,col-former2:col+latter2)))/0.3);
            if n1>=n2
                result1(row,col) = n1;
            else
                result2(row,col) = n2;
            end
        end
    end
end

%Show results
result3 = result1+result2; 
L_1vertical = result2;
L_2horizontal = result1; 
L = result3;

% Output steps lengths
ans_1vertical = sum(result2,"all")/xysize/512; %vertical
ans_2horizontal = sum(result1,"all")/xysize/512; %horizontal
ans_3_v_h = ans_1vertical/ans_2horizontal;
ans_4sum = sum(result3,"all")/xysize/512; %sum

% Gradient graphs superimpose step edge graphs
[I_grag,~] = imgradient(I_rescale,'sobel');
I_grag = rescale(I_grag);
I_grag = imgaussfilt(I_grag,2);
imshow(I_grag,[min(I_grag(:))*2 max(I_grag(:))*0.8])
hold on

% Divide orientation by Hough method
Angle_V = -45.1:0.01:45;
[H_V,T_V,R_V] = hough(result3,'RhoResolution',1.5,'Theta',Angle_V);
P_V = houghpeaks(H_V,999999,'threshold',ceil(0.3*max(H_V(:))));
lines_V = houghlines(result3,T_V,R_V,P_V,...
    'FillGap',5,'MinLength',(512/xysize)*0.85);
% Plot horizontal steps
for k = 1:length(lines_V)
   current_angle = lines_V(k).theta;
   if (0<=abs(current_angle)) && (abs(current_angle)<=5)
       xy_90 = [lines_V(k).point1; lines_V(k).point2];
       plot(xy_90(:,1),xy_90(:,2),'LineWidth',2,'Color','red');
   end
end

% Divide orientation by Hough method
Angle_H = [-90:0.01:-45 , 45.1:0.01:89.5];
[H_H,T_H,R_H] = hough(result3,'RhoResolution',1.5,'Theta',Angle_H);
P_H = houghpeaks(H_H,999999,'threshold',ceil(0.03*max(H_H(:))));
lines_H = houghlines(result3,T_H,R_H,P_H,...
    'FillGap',5,'MinLength',(512/xysize)*0.85);
% Plot vertical steps
for k = 1:length(lines_H)
   current_angle = lines_H(k).theta;
   if (85<=abs(current_angle)) && (abs(current_angle)<=90)
       xy_0 = [lines_H(k).point1; lines_H(k).point2];
       plot(xy_0(:,1),xy_0(:,2),'LineWidth',2,'Color','green');
   end
 end

lines = [lines_V,lines_H];
P = gcf;
exportgraphics(P,"f"+num2str(h)+".png")
hold off

% "Color" the edge extraction results.
imshow(zeros(512,512))
hold on

for k = 1:length(lines)
   current_angle = lines(k).theta;
   if (85<=abs(current_angle)) && (abs(current_angle)<=90)
       xy_90 = [lines(k).point1; lines(k).point2];
       plot(xy_90(:,1),xy_90(:,2),'LineWidth',2,'Color','green');
   elseif (0<=abs(current_angle)) && (abs(current_angle)<=5)
       xy_0 = [lines(k).point1; lines(k).point2];
       plot(xy_0(:,1),xy_0(:,2),'LineWidth',2,'Color','red');
   end
   
end
hold off

f = gcf;
exportgraphics(f,'f.png')
F = imread("f.png");
F = imresize(F,[512,512]);
F_V = F(:,:,1);
F_V(F_V>0) = 1;
F_H = F(:,:,2);
F_H(F_H>0) = 1;
I_V = result3.*double(F_V);
I_H = result3.*double(F_H);

for row = 1:512
    for col = 1:512
        if I_H(row,col) == 1
            [a1,pos1] = findbound(I_H,row,col,1);
%             [a2,pos2] = findbound(E_new,row,col,2);
            [former1,latter1] = findnei(a1,pos1);
%             [former2,latter2] = findnei(a2,pos2);
            n1 = round((max(I_raw(row-former1:row+latter1,col))-min(I_raw(row-former1:row+latter1,col)))/0.3);
%             n2 = round((max(I_raw(row,col-former2:col+latter2))-min(I_raw(row,col-former2:col+latter2)))/0.3);
                I_H(row,col) = n1;
%               I_V(row,col) = n2;

        end
    end
end

for row = 1:512
    for col = 1:512
        if I_V(row,col) == 1
%             [a1,pos1] = findbound(E_new,row,col,1);
            [a2,pos2] = findbound(I_V,row,col,2);
%             [former1,latter1] = findnei(a1,pos1);
            [former2,latter2] = findnei(a2,pos2);
%             n1 = round((max(I_raw(row-former1:row+latter1,col))-min(I_raw(row-former1:row+latter1,col)))/0.3);
            n2 = round((max(I_raw(row,col-former2:col+latter2))-min(I_raw(row,col-former2:col+latter2)))/0.3);
%                 I_H(row,col) = n1;
              I_V(row,col) = n2;

        end
    end
end

ans_11vertical = sum(I_V,"all")/xysize(1)/512; 
ans_22horizontal = sum(I_H,"all")/xysize/512; 
ans_33_v_h = ans_11vertical/ans_22horizontal;
ans_44sum = sum(I_V+I_H,"all")/xysize/512; 
ans_5lossR = (ans_4sum-ans_44sum)/ans_4sum;

end

% Function to get the length of SPM picture
function [length] = get_pic_length(filename) 
fileID = fopen(filename);
for i = 1:6
    tline = fgetl(fileID);
end
s=strsplit(tline);
length = str2double(s{end});
fclose(fileID);
end

% Function to get the filename of SPM picture
function [current_filename] = get_file_name(filename) 
fileID = fopen(filename);
for i = 1:3
    tline = fgetl(fileID);
end
s=strsplit(tline,"\");
current_filename = (s(end));
fclose(fileID);
end

function [array,pos] = findbound(E,y,x,direction)
% direction1:vertical 2:horizon
l=100;
if direction == 1
    if y <= l
        ystart = 1;
        yend = y+l;
        array = E(ystart:yend,x);
        pos = y;
    elseif y > size(E,1)-l
        ystart = y-l;
        yend = size(E,1);
        array = E(ystart:yend,x);
        pos = l+1;
    else
        ystart = y-l;
        yend = y+l;
        array = E(ystart:yend,x);
        pos = l+1;
    end
else
    if x <= l
        xstart = 1;
        xend = x+l;
        array = E(y,xstart:xend);
        pos = x;
    elseif x > size(E,2)-l
        xstart = x-l;
        xend = size(E,2);
        array = E(y,xstart:xend);
        pos = l+1;
    else
        xstart = x-l;
        xend = x+l;
        array = E(y,xstart:xend);
        pos = l+1;
    end
end
end


function [former,latter] = findnei(array,pos)
if length(find(array)) == 1
    former = floor(length(array(1:pos-1))/2);
    latter = ceil(length(array(pos+1:length(array)))/2);
elseif isempty(find(array, 1)) == 1
    former = 0;
    latter = 0;
else
    p1 = find(array(1:pos-1));
    p2 = find(array(pos+1:length(array)));
    if isempty(p1) == 1
        former = floor(length(array(1:pos-1))/2);
    else
        former = floor((pos - p1(length(p1)) - 1)/2);
    end
    if isempty(p2) == 1
        latter = ceil(length(array(pos+1:length(array)))/2);
    else
        latter = ceil((p2(1)-1)/2);
    end
end
end

function [U,S] = weighted_mean(x,w)
U = sum(x.*w)/sum(w);
S = sqrt( sum( w.*(x-U).^2 )/sum(w) );
end
