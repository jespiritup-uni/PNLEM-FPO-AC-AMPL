# Declaracion de conjuntos
set BAR;										# Conjunto de barras
set RAM within 1..500 cross BAR cross BAR;		# Conjunto de ramos

# Parametros de las barras
param tipo{BAR};          #Tipo de barra (1=PQ o carga, 2=PV o generacion, 3=Slack o referencia)
param Vb{BAR};            #Magnitud de tensión inicial o de operacion en la barra en p.u.
param angi{BAR};          #Angulo de la tension inical en la barra en grados
param Pd{BAR};            #Potencia Activa demandada en la barra
param Qd{BAR};            #Potencia Reactiva demandada en la barra                          
param Pgb{BAR};           #Potencia Activa generada en la barra
param Qgb{BAR};           #Potencia Reactiva generada en la barra                         
param gshb{BAR};          #Conductancia shunt en p.u. conectada a la barra
param bshb{BAR};          #Susceptancia shunt en p.u. conectada a la barra
param wshm{BAR};          #Estado del equipo de compensacion reactiva 1 si esta conectado y 0 si no esta conectado                         
param Pgm{BAR};           #Potencia Activa minima generada en la barra
param PgM{BAR};           #Potencia Activa maxima generada en la barra
param Qgm{BAR};           #Potencia Reactiva minima generada en la barra
param QgM{BAR};           #Potencia Reactiva maxima generada en la barra
param c2{BAR};            #coeficiente de la variable cuadratica de la funcion costo del generador
param c1{BAR};            #coeficiente de la variable lineal de la funcion costo del generador
param c0{BAR};            #coeficiente independiente del costo de generacion

# Parametros de los ramos
param r{RAM};             #Resistencia de la linea en p.u.
param x{RAM};             #Reactancia de la linea en p.u.
param bsh{RAM};           #Susceptancia shunt de la linea en p.u.
param fi{RAM};            #Angulo de desfase de transformadores cambiadores de fase en grados
param angmin{RAM};        #Angulo minimo permitido de desfase angular de la linea
param angmax{RAM};        #Angulo maximo permitido de desfase angular de la linea
param a{RAM};             #Relacion de transformacion del transformador
param SMAX{RAM};          #Capacidad maxima de flujo de potencia aparente de la linea de transmision
param ntap{RAM};          #Numero de posiciones discretos del Tap en el transformador
param rMtap{RAM};         #Relacion maxima de transformacion del transformador
param rmtap{RAM};         #Relacion minima de transformacion del transformador
                          
# parametros generales
param Sbase;              #Potencia base del sistema
param Vnom;               #Tension nominal p.u. del sistema
param Vmin;               #Tension nominal minimo en p.u. de operacion permitido para las barras
param Vmax;               #Tension nominal maximo en p.u. de operacion permitido para las barras

# Parametros adicionales para el FPO_CP
param g{RAM};             #Conductancia de la linea de transmision en p.u.
param b{RAM};             #Susceptancia de la linea de transmision en p.u.
param reg{RAM};           #Regulacion del tap del transformador en la linea

#Declaracion de variables de ramos
var nt{RAM} integer;      #Tap inicial del tranformador en la linea
var Ppara{RAM};           #Potencia activa de salida de la linea de transmision
var Qpara{RAM};           #Potencia reactiva de salida de la linea de transmision
var Pde{RAM};			  #Potencia activa de entrada de la linea de transmision
var Qde{RAM};             #Potencia reactiva de entrada de la linea de transmision

# Declaracion de las variables en las barras
var V{BAR};               #Tension de la barra
var th{BAR};              #Angulo de la tension en la barra
var Pger{BAR};            #Potencia activa generada en la barra
var Qger{BAR};            #Potencia reactiva generada en la barra
var gP{BAR} >= 0;         #Potencia activa en la barra
var gQ{BAR} >= 0;         #Potencia reactiva en la barra
var w{BAR} binary;        #Estado del Equipo de compensacion reactiva en la barra

#################################################
# Definicion de la funcion objetivo

minimize custo_generacion:
	sum {i in BAR: tipo[i] <> 0}(c2[i]*(Pger[i]*Sbase)^2 + c1[i]*Pger[i]*Sbase + c0[i]);

##################################################
# Definicion de las restriciones

subject to balance_potencia_activa_ex {i in BAR}:
	Pger[i] - Pd[i] - sum {(k,j,i) in RAM}(Ppara[k,j,i]) - sum {(k,i,j) in RAM}(Pde[k,i,j]) - V[i]^2*gshb[i] = 0;

subject to balance_potencia_reactiva_ex {i in BAR}:
	Qger[i] - Qd[i] - sum {(k,j,i) in RAM}(Qpara[k,j,i]) - sum {(k,i,j) in RAM}(Qde[k,i,j]) + V[i]^2*bshb[i]*w[i] = 0;

subject to flujo_potencia_activa_de {(k,i,j) in RAM}:
	Pde[k,i,j] = g[k,i,j]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]));

subject to flujo_potencia_activa_para {(k,i,j) in RAM}:
	Ppara[k,i,j] = g[k,i,j]*V[j]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]));

subject to flujo_potencia_reactiva_de {(k,i,j) in RAM}:
	Qde[k,i,j] = -(b[k,i,j] + bsh[k,i,j])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]));

subject to flujo_potencia_reactiva_para {(k,i,j) in RAM}:
	Qpara[k,i,j] = -(b[k,i,j] + bsh[k,i,j])*V[j]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]));

######################################################
# Definicion de los limites

# limite de tension
subject to limite_tension {i in BAR}:
	Vmin <= V[i] <= Vmax;

# limite de nt
subject to limite_nt {(k,i,j) in RAM}:
	- ntap[k,i,j] <= nt[k,i,j] <= ntap[k,i,j];

# limite do w
subject to limite_w {i in BAR}:
	0 <= w[i] <= wshm[i];

# Limite de la potencia P e Q
subject to limite_potencia_generada_activa {i in BAR}:
	Pgm[i] <= Pger[i] <= PgM[i];

subject to limite_potencia_generada_reactiva {i in BAR}:
	Qgm[i] <= Qger[i] <= QgM[i];
	
# cuadrado de flujo de potencia P e Q
subject to cuadrado_flujo_potencia_de {(k,i,j) in RAM}:
	Pde[k,i,j]^2 + Qde[k,i,j]^2 <= SMAX[k,i,j]^2;

subject to cuadrado_flujo_potencia_para {(k,i,j) in RAM}:
	Ppara[k,i,j]^2 + Qpara[k,i,j]^2 <= SMAX[k,i,j]^2;
