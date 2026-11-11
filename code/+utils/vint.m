function y = vint(F, w)
%VINT Weighted velocity integral. F is Nv-by-Nx or Nv-by-1, w is Nv-by-1.
w = w(:);
y = w.' * F;
end
