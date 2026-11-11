function y = phi_response(z, p)
switch lower(p.phi_type)
    case 'sign'
        y = -sign(z);
    case 'exp_sign'
        % Response alpha is independent of p.alpha (chemical degradation).
        assert(isfield(p,'phi_alpha') && isscalar(p.phi_alpha) ...
            && isfinite(p.phi_alpha) && p.phi_alpha>=0, ...
            'exp_sign requires a finite nonnegative p.phi_alpha.');
        y = -exp(-p.phi_alpha*abs(z)).*sign(z);
    case 'linear'
        y = -z;
    case 'tanh'
        y = -tanh(z/max(p.phi_delta, eps));
    otherwise
        error('Unknown response type: %s',p.phi_type);
end
end
