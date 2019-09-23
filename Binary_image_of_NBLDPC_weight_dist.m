
%Optimization of NB-LDPC row codes
%Use part of WeightDistribution.m from % Copyright (c) 2006. Robert Morelos-Zaragoza. All rights reserved.
%Thank You Robert for your amazing practical related book. Your Source Code support 
%several generation of Channel Coding Scientist!

%Copyright(c) 2012, USAtyuk Vasiliy 
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met :
%*Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and / or other materials provided with the distribution.
%* Neither the name of the <organization> nor the
%names of its contributors may be used to endorse or promote products
%derived from this software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%logq=4;
%u=[1 1 0 0];%primitive polynomial coefficient vector starting p0 to p_(field-1)

logq=8;
u=[1 0 1 1 1 0 0 0];%primitive polynomial coefficient vector starting p0 to p_(field-1)
[row, column  ]=size(H);

powers=0:2^logq-1;
%powers=randi([2^logq-1],1,6);
M=[[zeros(logq-1,1),eye(logq-1)]];
M=vertcat(M,u);%generate companion matrix
Binary_image(logq, logq, size(powers,2))=zeros();
for i=1:size(powers,2)
Binary_image(:,:,i)=M;
end
for i=1:size(powers,2)
if powers(i)==0
Binary_image(:,:,i)=eye(logq);
end
if powers(i)~=0
for j=1:powers(i)-1
  Binary_image(:,:,i)=mod(Binary_image(:,:,i)*M,2);
end
end
end
for i=1:size(powers,2)
Binary_image(:,:,i)=Binary_image(:,:,i)';
end


HBinImage=[];
for i=1:row
    HBinImagerow=[];
    for j=1:column
        if(H(i,j)==0)
            HBinImagerow=horzcat(HBinImagerow,Binary_image(:,:,1));
        end
        if(H(i,j)==-1)
            HBinImagerow=horzcat(HBinImagerow,zeros(logq,logq));
        end
        if(H(i,j)~=0)&&(H(i,j)~=-1)
            HBinImagerow=horzcat(HBinImagerow,Binary_image(:,:,H(i,j)+1));
        end
    end
    HBinImage=vertcat(HBinImage,HBinImagerow);
end

HBinImage=sparse(HBinImage);

[h, g] = h2g(HBinImage,2) ;
sparse(mod(g*h',2))

%estimate weight distribution



%n=10;               % Code length
%n=8;               % Code length
%k=4;                % Code dimension
% Generator matrix, k=4, n=10
% G = [ 1 1 0 1 0 0 0 0 0 0
%       0 0 1 1 0 1 0 0 0 0
%       0 0 0 0 1 1 0 1 0 0
%       0 0 0 0 0 0 1 1 0 1 ];
% G = [ 1 1 0 1 0 0 0 0 
%       0 0 1 1 0 1 0 0 
%       0 0 0 0 1 1 0 1 
%       0 0 0 0 0 0 1 1 ];
%G = [ 1 1 0 1 0 0 0 0 
%      0 0 1 1 0 1 0 0 
%      0 0 0 0 1 1 0 1 
%      0 1 0 0 0 0 1 1 ];
%G = [[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0;0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]];    

n=size(g,2); %Code length
k=size(g,1); % Code dimension

% Weight distribution W(C):
W = zeros(1,n+1);
W(1)=1;

% Generate all combinations of message vectors and their weigths
k2=2^k-1;
for i=1:k2
    u = dec2bin(i,k);
    v = mod(u*g,2);
    W(sum(v)+1) = W(sum(v)+1)+1;
end

% Print the weight distribution
fprintf('W(C)={1');
for i=2:n+1
    fprintf(',%d',W(i))
end
fprintf('}\n');

% Compute the union bound

EbNo = 0:0.5:9.5;
EbNoratio = 10.^(EbNo/10);
bound = 0;
for i=1:n
    bound = bound + i*W(i+1)/n * erfc(sqrt(2*i*(k/n)*EbNoratio));
end

% semilogy(EbNo,bound,'-o')
semilogy(EbNo,bound,'-s')
axis tight
grid on
xlabel('E_b/N_0 (dB)')
ylabel('P_b')
hold on

% semilogy(EbNo, Q(sqrt(2*EbNoratio)))