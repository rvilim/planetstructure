def K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma):
	if T is False:
		return K_T(rho, rho_0, T, K_0, Kp_0)*(1+alpha*gamma*300.0)

	return K_T(rho, rho_0, T, K_0, Kp_0)*(1+alpha*gamma*T)

def K_T(rho, rho_0, T, K_0, Kp_0):
	if T is False:
		return K_T_rho300(rho, rho_0, K_0, Kp_0)
	else:
		return K_T_rho300(rho, rho_0, K_0, Kp_0)+DeltaK_th(rho, T)

def K_T_rho300(rho, rho_0, K_0, Kp_0):
	x=rho/rho_0
	# return 3*K_0*(pow(x,2.0/3)-pow(x,1.0/3))*math.exp(1.5*(Kp_0-1)*(1-pow(x,-1.0/3)))
	return (K_0/2.0)*(7*pow(x,7.0/3.0)-5*pow(x,5.0/3.0)*(1+.75*(Kp_0-4)*(pow(x,2.0/3)-1))+1.5*(pow(x,3.0)-pow(x,7.0/3))*(Kp_0-4.0))

def DeltaK_th(rho, T):
	return 3*n*R*gamma*rho*(f(T)-f(T_0))

def f_integrand(xi):
	return pow(xi,3)/(math.exp(xi)-1)+3*theta*gamma/(math.exp(theta/T)-1)

def f(T):
	theta=theta0*math.exp((gamma0-gamma)/q)

	return (1-q-3*gamma)*(pow(T,4)/pow(theta,3))*integrate.quad(f_integrand,0,theta/T)