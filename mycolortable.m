function CT = mycolortable(varargin)
% CT = MYCOLORTABLE(N)
%   Return a color table
%
%  N optionally specifies the length of the table
%  (default 256)
%

N = 256;
if numel(varargin)>0
	thisarg = varargin{1};
	if isnumeric(thisarg) && isscalar(thisarg)
		N = thisarg;
	end
end

CT0 = [0.1405 0.00719 0.2242;0.2134 0.02435 0.3071;0.2648 0.02239 0.3479;0.2939 0.02786 0.3459;0.3375 0.02245 0.3691;0.3897 0.0223 0.4038;0.4351 0.02615 0.4378;0.456 0.01441 0.4381;0.4954 0.01003 0.449;0.5511 0.009804 0.4672;0.6116 0.008296 0.464;0.6801 0.001307 0.4526;0.7271 0.00228 0.4565;0.7818 0.006536 0.4948;0.828 0.009907 0.5359;0.8676 0.0143 0.5839;0.8606 0.01917 0.6722;0.8253 0.02237 0.7799;0.818 0.02273 0.8655;0.7843 0.01825 0.9707;0.7034 0.02612 0.9948;0.6356 0.02984 0.975;0.5591 0.0158 0.986;0.4729 0.01941 0.9852;0.3943 0.01963 0.9843;0.2488 0.0183 0.9863;0.07647 0.02211 0.9914;0.008018 0 0.9843;0.01111 0 0.969;0.01767 0.01238 0.9275;0.01273 0.01356 0.8648;0.03213 0.01246 0.8417;0.02846 0.04258 0.801;0.0274 0.1747 0.8174;0.02555 0.2225 0.8204;0.03333 0.2876 0.817;0.02546 0.3374 0.8016;0.0236 0.3879 0.8185;0.01625 0.432 0.8365;0.006413 0.4712 0.8649;0.007416 0.4887 0.8925;0.0007977 0.5012 0.9207;0.0004992 0.5309 0.9578;0 0.5756 0.9916;0.003567 0.6182 0.997;0.00719 0.6767 0.9928;0.02484 0.7205 0.9989;0.03152 0.7562 0.993;0.01508 0.822 0.9863;0.01747 0.8936 0.9953;0.02296 0.9529 0.9891;0.024 0.9852 0.9279;0.02026 0.9961 0.8719;0.01496 0.9886 0.7667;0.03049 0.9892 0.6001;0.02087 0.9863 0.4098;0.01789 0.9974 0.07776;0.01725 0.9606 0.0003963;0.01863 0.9308 0.01073;0.001961 0.9131 0;0.004236 0.8942 0;0.005677 0.8781 0.001307;0.005311 0.8535 0.0006948;0.01288 0.833 0.01111;0.01661 0.8162 0.01895;0.01176 0.7963 0.01565;0.02207 0.786 0.02584;0.02157 0.7644 0.02651;0.06474 0.7566 0.03014;0.1571 0.7701 0.01745;0.1626 0.7964 0.02131;0.1612 0.8166 0.008401;0.2329 0.8254 0.01182;0.2985 0.8546 0;0.3777 0.8741 0.01526;0.474 0.8925 0.006012;0.5589 0.9105 0.0002779;0.6113 0.9272 0.01176;0.6918 0.9392 0.01242;0.7733 0.9464 0.02611;0.8095 0.9646 0.02408;0.87 0.9698 0.03896;0.9162 0.9704 0.0406;0.9458 0.9832 0.02098;0.9722 0.9876 0.03364;0.9933 0.9877 0.03591;0.9966 0.9901 0.03054;0.9943 0.9878 0.0268;0.9928 0.9864 0.02418;0.9907 0.9778 0.01195;0.9969 0.9605 0.01176;0.9987 0.9372 0.01176;0.9908 0.9186 0.004179;0.987 0.8869 0.005882;0.9922 0.8495 0.007843;0.9884 0.8197 0.01942;0.9843 0.7858 0.02226;0.9892 0.7445 0.02651;0.9838 0.7106 0.01977;0.9882 0.677 0.01548;0.9961 0.6406 0.02876;0.9924 0.6112 0.02539;0.9926 0.583 0.03174;0.9967 0.5562 0.02285;0.999 0.538 0.01752;0.9943 0.523 0.00915;0.9987 0.5125 0.009305;0.9975 0.497 0;0.9974 0.4951 0;0.9984 0.49 0.00244;1 0.4782 0.001307;0.9978 0.4522 0.003278;0.9975 0.434 0.01569;0.9972 0.4088 0.01825;0.998 0.3756 0.02078;0.9958 0.3408 0.01079;0.9882 0.3182 0.003;0.9839 0.2809 0.01128;0.9916 0.2318 0.01659;0.9952 0.1948 0.01917;0.9895 0.1604 0.01373;0.9863 0.1382 0.01373;0.9863 0.09212 0.002857;0.9938 0.05828 0.003598;0.9942 0.04339 0.01283;0.9941 0.02184 0.01438;0.9882 0.02519 0.01642;0.985 0.02092 0.03268];

na = size(CT0,1);
CT = interp1(linspace(0,1,na),CT0,linspace(0,1,N));
