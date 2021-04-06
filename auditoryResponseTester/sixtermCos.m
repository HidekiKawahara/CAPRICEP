function output = sixtermCos(n_length)
%  output = sixtermCos(n_length)
%  generates the six-term cosine series with very low side lobes and steep
%  decay rate
%
% Argument
%   n_length  : length of the generated shape (samples)
%
% Output
%   output    : generated cosine series
%
% Reference
%   [1] Kawahara, H., Sakakibara, K., Morise, M., Banno, H., Toda, T., 
%   Irino, T. (2017) A New Cosine Series Antialiasing Function and its 
%   Application to Aliasing-Free Glottal Source Models for Speech and 
%   Singing Synthesis. Proc. Interspeech 2017, 1358-1362, 
%   DOI: 10.21437/Interspeech.2017-15.
%
% Coded by Hideki Kawahara 18 February 2018: creation date
% Attribution 4.0 International (CC BY 4.0)
% https://creativecommons.org/

base_index = ((0:n_length - 1)' - n_length / 2 + 0.5) / n_length;
BB = [0.2624710164;0.4265335164;0.2250165621;0.0726831633;0.0125124215;0.0007833203];
output = cos(base_index * [0 1 2 3 4 5] * 2 * pi) * BB;
end