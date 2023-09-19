function [w] = weightupdate(pD,absW,beta_FA,d,prev,estv,preV,estV,hat_X,hat_R,S,preAlpha,estAlpha,preBeta,estBeta)
% B.Reference:[1] K. Granstrom, A. Natale, P. Braca, G. Ludeno, and F. Serafino, 
% "Gamma Gaussian Inverse Wishart Probability Hypothesis Density for Extended Target Tracking Using X-Band Marine Radar Data," 
% IEEE Transactions on Geoscience and Remote Sensing.


part1 = absW*d/2*log(2*pi) + absW/2*log(2)-d/2*log(absW);   %(1st part of 38e)
part2 = prev/2*log(det(preV)) - estv/2*log(det(estV));
part3 = gammaln(estv/2) - gammaln(prev/2);
part4 = absW/2*log(det(hat_X)) - (absW-1)/2*log(det(hat_R)) - 1/2*log(det(S));
part5 = preAlpha*log(gamma(estAlpha)*preBeta) - estAlpha*log(gamma(preAlpha)*estBeta);
logLikelihood = part1 + part2 + part3 + part4 + part5;
w = log(pD) + logLikelihood - absW*log(beta_FA);





end