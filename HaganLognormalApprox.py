# -*- coding: utf-8 -*-
import math

def haganLogNormalApprox(y, expiry, F_0, alpha_0, beta, nu, rho):
    """Function which returns the Black implied volatility, computed using the Hagan et al. 
    lognormal approximation"""
    
    one_beta = 1.0 - beta
    one_betasqr = one_beta * one_beta
    if F_0!=y:
        fK = F_0 *y
        fK_beta = math.pow(fK, one_beta/ 2.0)
        log_fK = math.log(F_0/y)
        z = nu/alpha_0  *fK_beta * log_fK
        x = math.log((math.sqrt(1.0-2.0*rho*z +z*z)+ z -rho)/(1-rho))
        sigma_1 = (alpha_0 / fK_beta / (1.0 + one_betasqr / 24.0 *log_fK *log_fK + math.pow(one_beta *log_fK,4)/1920.0)*(z/x))
        sigma_exp = (one_betasqr / 24.0 * alpha_0*alpha_0/fK_beta /fK_beta + 0.25 *rho *beta*nu*alpha_0/fK_beta + (2.0 - 3.0 * rho * rho)/ 24.0 *nu *nu) 
        sigma =sigma_1 * (1.0 + sigma_exp*expiry)
    else:
        f_beta = math.pow(F_0,one_beta)
        f_two_beta = math.pow(F_0,(2.0 - 2.0 *beta))
        sigma = ((alpha_0/f_beta)*(1.0 + ((one_betasqr/24.0)*(alpha_0 *alpha_0/f_two_beta)+(0.25*rho*beta*nu*alpha_0/f_beta)+ (2.0 -3.0 *rho*rho)/24.0*nu*nu)*expiry))
        
    return sigma