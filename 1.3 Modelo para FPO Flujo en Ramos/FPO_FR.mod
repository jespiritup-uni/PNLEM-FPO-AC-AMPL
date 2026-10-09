# Declaracion de los conjuntos
set BAR;										# Conjunto de barras
set RAM within 1..500 cross BAR cross BAR;		# Conjunto de ramos

# Parametros de las barras
param tipo{BAR};     #Tipo de barra (1=PQ o carga, 2=PV o generacion, 3=Slack o referencia)
param Vb{BAR};       #Magnitud de tensión inicial o de operacion en la barra en p.u.
param angi{BAR};     #Angulo de la tension inical en la barra en grados
param Pd{BAR};       #Potencia Activa demandada en la barra
param Qd{BAR};       #Potencia Reactiva demandada en la barra				     
param Pgb{BAR};      #Potencia Activa generada en la barra
param Qgb{BAR};      #Potencia Reactiva generada en la barra				     
param gshb{BAR};     #Conductancia shunt en p.u. conectada a la barra
param bshb{BAR};     #Susceptancia shunt en p.u. conectada a la barra
param wshm{BAR};     #Estado del equipo de compensacion reactiva 1 si esta conectado y 0 si no esta conectado				     
param Pgm{BAR};      #Potencia Activa minima generada en la barra
param PgM{BAR};      #Potencia Activa maxima generada en la barra
param Qgm{BAR};      #Potencia Reactiva minima generada en la barra
param QgM{BAR};      #Potencia Reactiva maxima generada en la barra
param c2{BAR};       #coeficiente de la variable cuadratica de la funcion costo del generador
param c1{BAR};       #coeficiente de la variable lineal de la funcion costo del generador
param c0{BAR};       #coeficiente independiente del costo de generacion

# Parametros de los ramos
param r{RAM};        #Resistencia de la linea en p.u.
param x{RAM};		 #Reactancia de la linea en p.u.
param bsh{RAM};      #Susceptancia shunt de la linea en p.u.
param fi{RAM};       #Angulo de desfase de transformadores cambiadores de fase en grados
param angmin{RAM};   #Angulo minimo permitido de desfase angular de la linea
param angmax{RAM};   #Angulo maximo permitido de desfase angular de la linea
param a{RAM};        #Relacion de transformacion del transformador
param SMAX{RAM};	 #Capacidad maxima de flujo de potencia aparente de la linea de transmision
param ntap{RAM};	 #Numero de posiciones discretos del Tap en el trasformador
param rMtap{RAM};	 #Relacion maxima de transformacion del transformador
param rmtap{RAM};	 #Relacion minima de transformacion del transformador

# parametros generales
param Sbase;         #Potencia base del sistema
param Vnom;			 #Tension nominal p.u. del sistema
param Vmin;			 #Tension nominal minimo en p.u. de operacion permitido para las barras
param Vmax;			 #Tension nominal maximo en p.u. de operacion permitido para las barras

# Parametros adicionales para el FPO_FR
param z2{RAM};       #Impedancia elevado al cuadrado de la linea en p.u.
param reg{RAM};      #Regulacion del tap del transformador en la linea

# Declaracion de las variables en las barras
var Vqdr{BAR};       #Tension elevado al cuadrado de la barra
var theta{BAR};      #Angulo de la tension de la barra				     
var Pger{BAR};	     #Potencia activa generada en la barra
var Qger{BAR};	     #Potencia reactiva generada en la barra
var w{BAR} binary;   #Estado del Equipo de compensacion reactiva en la barra

# Declaracion de las variables en los ramos
var Iqdr{RAM} >= 0;  #Corriente elevado al cuadrado de la linea
var nt{RAM} integer; #Tap inicial del transformador en la linea
var P{RAM};          #Potencia activa de la linea
var Q{RAM};			 #Potencia reactiva de la linea

#################################################
# Definicion de la funcion objetivo

minimize custo_generacion:
	sum {i in BAR: tipo[i] <> 0}(c2[i]*(Sbase*Pger[i])^2 + c1[i]*Sbase*Pger[i] + c0[i]);

##################################################
# Definicion de las restricciones

subject to balance_potencia_activa {i in BAR}:
	Pger[i] - Pd[i] + sum {(k,j,i) in RAM}(P[k,j,i]) - sum {(k,i,j) in RAM}(P[k,i,j] + r[k,i,j]*Iqdr[k,i,j]) - Vqdr[i] * gshb[i]= 0;

subject to balance_potencia_reactiva {i in BAR}:
	Qger[i] - Qd[i] + sum {(k,j,i) in RAM}(Q[k,j,i] + bsh[k,j,i]*Vqdr[i]) - sum {(k,i,j) in RAM}(Q[k,i,j] - bsh[k,i,j]*Vqdr[i] + x[k,i,j]*Iqdr[k,i,j]) + Vqdr[i]*bshb[i]*w[i] = 0;

subject to caida_de_tension_linea {(k,i,j) in RAM}:
	Vqdr[i]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2 - Vqdr[j] = 2*(r[k,i,j]*P[k,i,j] + x[k,i,j]*Q[k,i,j]) + z2[k,i,j]*Iqdr[k,i,j];

subject to abertura_angular_linea_exacta {(k,i,j) in RAM}:
	sqrt(Vqdr[i])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*sqrt(Vqdr[j])*sin(theta[i] - theta[j] + fi[k,i,j]) = x[k,i,j]*P[k,i,j] - r[k,i,j]*Q[k,i,j];

subject to calculo_magnitud_corriente_exacta {(k,i,j) in RAM}:
	Vqdr[j]*Iqdr[k,i,j] = P[k,i,j]^2 + Q[k,i,j]^2;

######################################################
# Definicion de los limites

# limite de tension
subject to limite_tension {i in BAR}:
	Vmin^2 <= Vqdr[i] <= Vmax^2;

# limite de nt
subject to limite_nt {(k,i,j) in RAM}:
	-ntap[k,i,j] <= nt[k,i,j] <= ntap[k,i,j];

# limite do w
subject to limite_w {i in BAR}:
	0 <= w[i] <= wshm[i];

# Limite do potencia P e Q
subject to limite_potencia_generada_activa {i in BAR}:
	Pgm[i] <= Pger[i] <= PgM[i];

subject to limite_potencia_generada_reactiva {i in BAR}:
	Qgm[i] <= Qger[i] <= QgM[i];

subject to cuadrado_flujo_potencia_1 {(k,i,j) in RAM}:
	P[k,i,j]^2 + (Q[k,i,j] + bsh[k,i,j]*Vqdr[i])^2 <= SMAX[k,i,j]^2;

subject to cuadrado_flujo_potencia_2 {(k,i,j) in RAM}:
	(P[k,i,j] + r[k,i,j]*Iqdr[k,i,j])^2 + (Q[k,i,j] - bsh[k,i,j]*Vqdr[i] + x[k,i,j]*Iqdr[k,i,j])^2 <= SMAX[k,i,j]^2;
