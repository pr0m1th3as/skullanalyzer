% Copyright (C) 2019 Andreas Bertsatos <andreas.bertsatos@gmail.com>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

function plot_features(varargin)
	% -*- texinfo -*-
  % @deftypefn  { plot_features } 
  % @deftypefnx  {plot_features(@var{filename})} 
  %
  % This function reads the available csv files, which contain the Nasion-Bregma 2D polyline
  % and other HMIs extracted with the 'skullanalyzer' program, and plots them in separate
	% figures. It requires the base filename as an input argument or it prompts the user in its
	% absense. The function assumes that all csv files use the naming convention defined by
	% 'skullanalyzer' and as such, the base filename should be the filename of the analyzed 3D
	% model without its extension (.obj). Furthermore, each plot is saved in a .png image following
	% an appropriate name convention.
	%
  % The present function is complementary to the 'skullanalyzer' program which can be found at
  % @url{https://github.com/pr0m1th3as/skullanalyzer}.
	% @end deftypefn
	
  % check for input argument or ask user for filename
  if(nargin == 1)
    filename = varargin{1};
  else
    filename = inputdlg("Enter Mesh name");
  endif
  extension = ["*.csv"];
  filenames = ls(char(strcat(filename, extension)));
  name = strcat(filename,".obj");
  % iterate through the list of available files and plot them in separate figures
  for i=1:length(filenames(:,1))
    % for each filename load data in particular variable
    % Nasion Bregma segment
    if(index(filenames(i,:), strcat(filename,"_NasionBregmaSegment")) == 1)
      NasionBregmaSegment = csvread(strcat(filenames(i,:)));
    endif
		% occipital protuberance Height Map Image
		if(index(filenames(i,:), strcat(filename,"_occipitalHeightMapImage.csv")) == 1)
			occipitalHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
		% supraorbital ridge Height Map Image
		if(index(filenames(i,:), strcat(filename,"_supraorbitalHeightMapImage.csv")) == 1)
			supraorbitalHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
		% left mastoid process Lateral Height Map Image
		if(index(filenames(i,:), strcat(filename,"_mastoidLeftLateralHeightMapImage.csv")) == 1)
			mastoidLeftLateralHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
		% left mastoid process Inferior Height Map Image
		if(index(filenames(i,:), strcat(filename,"_mastoidLeftInferiorHeightMapImage.csv")) == 1)
			mastoidLeftInferiorHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
		% right mastoid process Lateral Height Map Image
		if(index(filenames(i,:), strcat(filename,"_mastoidRightLateralHeightMapImage.csv")) == 1)
			mastoidRightLateralHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
		% right mastoid process Inferior Height Map Image
		if(index(filenames(i,:), strcat(filename,"_mastoidRightInferiorHeightMapImage.csv")) == 1)
			mastoidRightInferiorHeightMapImage = csvread(strcat(filenames(i,:)));
		endif
  endfor
  % plot figures to specified order
  % Nasion Bregma segment
  if(exist("NasionBregmaSegment"))
    h1 = figure;
    plot(NasionBregmaSegment(:,1),NasionBregmaSegment(:,2), 'linestyle', '-', 'marker', '.');
    axis("equal");
    text = ["Nasion Bregma segment of " name "\nRight hand side view"];
    title(text);
		print(h1, "NasionBregma.png");
  endif
	% occipital protuberance Height Map Image
	if(exist("occipitalHeightMapImage"))
		h2 = figure;
		imshow(occipitalHeightMapImage, [0, 255]);
		text = ["Occipital Protuberance height map of " name "\nPosterior view"];
    title(text);
		print(h2, "OccipitalProtuberance.png");
	endif
	% supraorbital ridge Height Map Image
	if(exist("supraorbitalHeightMapImage"))
		h3 = figure;
		imshow(supraorbitalHeightMapImage, [0, 255]);
		text = ["Supraorbital Ridge height map of " name "\nAnterior view"];
    title(text);
		print(h3, "SupraorbitalRidge.png");
	endif
	% left mastoid process lateral Height Map Image
	if(exist("mastoidLeftLateralHeightMapImage"))
		h4 = figure;
		imshow(mastoidLeftLateralHeightMapImage, [0, 255]);
		text = ["Left Mastoid Process height map of " name "\nLateral view"];
    title(text);
		print(h4, "mastoidLeftLateral.png");
	endif
	% left mastoid process inferior Height Map Image
	if(exist("mastoidLeftInferiorHeightMapImage"))
		h5 = figure;
		imshow(mastoidLeftInferiorHeightMapImage, [0, 255]);
		text = ["Left Mastoid Process height map of " name "\nInferior view"];
    title(text);
		print(h5, "mastoidLeftInferior.png");
	endif
	% right mastoid process lateral Height Map Image
	if(exist("mastoidRightLateralHeightMapImage"))
		h6 = figure;
		imshow(mastoidRightLateralHeightMapImage, [0, 255]);
		text = ["Right Mastoid Process height map of " name "\nLateral view"];
    title(text);
		print(h6, "mastoidRightLateral.png");
	endif
	% right mastoid process inferior Height Map Image
	if(exist("mastoidRightLateralHeightMapImage"))
		h7 = figure;
		imshow(mastoidRightInferiorHeightMapImage, [0, 255]);
		text = ["Right Mastoid Process height map of " name "\nInferior view"];
    title(text);
		print(h7, "mastoidRightInferior.png");
	endif
endfunction
