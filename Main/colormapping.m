function[c1,n_div,div] = colormapping(Phi,n)
% colors for phi
div_max=90;
div=-div_max:(5):div_max;
n_div=length(div)-1;
color=jet(n_div);

for i = 1:n
    % color mapping for Phi angle
    c1(i,:)=color(find(div>Phi(i),1)-1,:);
end

end