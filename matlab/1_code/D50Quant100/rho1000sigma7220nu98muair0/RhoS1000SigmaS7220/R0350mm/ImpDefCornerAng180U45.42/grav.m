function grav=grav(Gamma,t,g,wzero,thetaZero)
grav = g*(1-Gamma*cos(wzero*t+thetaZero));