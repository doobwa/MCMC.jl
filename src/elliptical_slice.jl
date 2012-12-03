
# http://homepages.inf.ed.ac.uk/imurray2/pub/10ess/elliptical_slice.m
# https://github.com/jobovy/bovy_mcmc/blob/master/bovy_mcmc/elliptical_slice.py


function elliptical_slice(xx, prior::Array{Float64,}, log_like_fn, cur_log_like)
    D = length(xx);
    cur_log_like = log_like_fn(xx);
    angle_range = 0;
    hh = log(rand()) + cur_log_like;

    # Set up a bracket of angles and pick a first proposal.
    # "phi = (theta'-theta)" is a change in angle.
    if angle_range <= 0
        # Bracket whole ellipse with both edges at first proposed point
        phi = rand() * 2*pi;
        phi_min = phi - 2*pi;
        phi_max = phi;
    else
        # Randomly center bracket on current point
        phi_min = -angle_range*rand();
        phi_max = phi_min + angle_range;
        phi = rand()*(phi_max - phi_min) + phi_min;
    end

    xx_prop = xx  # TODO: Need this due to scoping rules?

    # Slice sampling loop
    while true
        # Compute xx for proposed angle difference and check if it's on the slice
        xx_prop = xx*cos(phi) + nu*sin(phi);
        cur_log_like = log_like_fn(xx_prop);
        if cur_log_like > hh
            # New point is on slice, ** EXIT LOOP **
            break;
        end
        # Shrink slice to rejected point
        if phi > 0
            phi_max = phi;
        elseif phi < 0
            phi_min = phi;
        else
            error("BUG DETECTED: Shrunk to current position and still not acceptable.");
        end
        # Propose new angle difference
        phi = rand() * (phi_max - phi_min) + phi_min;
    end
    return xx_prop, cur_log_like
end

# User specified Cholesky of prior covariance:
function elliptical_slice(xx, Sigma::Array{Float64,2}, llk, cur_llk)
    nu = reshape(Sigma'*randn(D, 1), size(xx));
    elliptical_slice(xx, nu, llk, cur_llk)
end
