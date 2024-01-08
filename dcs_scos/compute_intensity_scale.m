function [intDist,intScale]=compute_intensity_scale(dataDir,ref_distance,ref_CR)
%% Computes the scaling factor of the intensity based on the Monte Carlo simulation
%  and a reference count rate found from a DCS measurement

% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% dataDir: Directory where the subject model and results are saved. Works
% with defined behaviour when there are only one of each file type of
% interest is in the folder.

% Outputs
% intDist: Source-detector separations given in millimeters from the source
% intScale: Scale factor for the count rate for a given distance from the
% source

arguments
dataDir char
ref_distance (1,1) double
ref_CR (1,1) double
end

%% Find all the necessary files
files=dir(dataDir);

% json input file
selFile=contains({files.name},'.json');
jsonFile=files(selFile);

% diffuse reflectance file
selFile=contains({files.name},'_dref.dat');
drefFile=files(selFile);

% Find the node file
selFile=and(contains({files.name},'node_'),contains({files.name},'.dat'));
nodeFile=files(selFile);

% Find the face file
selFile=and(contains({files.name},'face_'),contains({files.name},'.dat'));
faceFile=files(selFile);

%% Load the inputs

% Load the json file to get the source position
jsonFID=fopen([dataDir filesep jsonFile.name],'r');
jsonText=fscanf(jsonFID,'%s',[1,inf]);
jsonStruct=jsondecode(jsonText);
fclose(jsonFID);
srcPos=jsonStruct.Optode.Source.Pos';
detPos=[jsonStruct.Optode.Detector.Pos]';

% Load the diffuse reflectance file
drefFID=fopen([dataDir filesep drefFile.name],'r');
dref=fscanf(drefFID,'%d %e',[2,inf]);
fclose(drefFID);

% Load the node and face positions
node=readmmcnode([dataDir filesep nodeFile.name]);
face=readmmcface([dataDir filesep faceFile.name]);
avgNode=mean(node,1);

%% Compute the diffuse reflectance as a function of source-detector separation

% Select the points that are non-zero and less than 6 cm from the source
selDref=dref(2,:)>0;

euDist=sqrt((srcPos(1)-mean(reshape(node(face',1),[3,size(face,1)]),1)).^2+...
    (srcPos(2)-mean(reshape(node(face',2),[3,size(face,1)]),1)).^2+...
    (srcPos(3)-mean(reshape(node(face',3),[3,size(face,1)]),1)).^2);
selDist=euDist<=60;

selIndex=and(selDref,selDist);

% Give the settings for the averaged curve
blurStd=1;
intDist=0:.1:60;
weight=zeros(size(intDist));
scale=zeros(size(intDist));

% Select the face and diffuse reflectance information from the full list
selFace=face(selIndex,:);
selDRef=dref(2,selIndex);

% Compute the surface distance between the source and the data point
faceD=zeros([1,sum(selIndex)]);
faceA=zeros([1,sum(selIndex)]);
for faceIter=1:size(selFace,1)
    euDist=sqrt((srcPos(1)-mean(node(selFace(faceIter,:),1))).^2+...
    (srcPos(2)-mean(node(selFace(faceIter,:),2))).^2+...
    (srcPos(3)-mean(node(selFace(faceIter,:),3))).^2);
    if euDist>=3
    [faceD(faceIter),~]=find_spherical_distance(node,face,srcPos,...
        [mean(node(selFace(faceIter,:),1)),mean(node(selFace(faceIter,:),2)),mean(node(selFace(faceIter,:),3))],...
        avgNode);
    else
        faceD(faceIter)=euDist;
    end
    vec1=[node(selFace(faceIter,2),1)-node(selFace(faceIter,1),1),...
        node(selFace(faceIter,2),2)-node(selFace(faceIter,1),2),...
        node(selFace(faceIter,2),3)-node(selFace(faceIter,1),3)];
    vec2=[node(selFace(faceIter,3),1)-node(selFace(faceIter,1),1),...
        node(selFace(faceIter,3),2)-node(selFace(faceIter,1),2),...
        node(selFace(faceIter,3),3)-node(selFace(faceIter,1),3)];
    faceA(faceIter)=sqrt(sum(cross(vec1,vec2).^2))/2;
    
    weight=weight+exp(-(faceD(faceIter)-intDist).^2/(2*blurStd^2));
    scale=scale+exp(-(faceD(faceIter)-intDist).^2/(2*blurStd^2))*log10(selDRef(faceIter)/faceA(faceIter));
end
intScale=10.^(scale./weight)./mean(10.^(scale./weight));


sds=zeros([1,size(detPos,1)]);
for sdsIter=1:length(sds)
    [sds(sdsIter),~]=find_spherical_distance(node,face,srcPos,detPos(sdsIter,:),avgNode);
end
intScale=interp1(intDist,intScale,sds)./interp1(intDist,intScale,ref_distance)*ref_CR;
intDist=sds;
end
