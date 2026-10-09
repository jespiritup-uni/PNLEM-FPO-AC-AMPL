# Declaracion de conjuntos
set BAR;										# Conjunto de barras
set RAM within 1..500 cross BAR cross BAR;		# Conjunto de ramos

# Parametros de las BARRAS
param tipo{BAR};				#Tipo de barra (1=PQ o carga, 2=PV o generacion, 3=Slack o referencia)
param Vb{BAR};                  #Magnitud de tensión inicial o de operacion en la barra en p.u.
param angi{BAR};                #Angulo de la tension inical en la barra en grados
param Pd{BAR};                  #Potencia Activa demandada en la barra
param Qd{BAR};                  #Potencia Reactiva demandada en la barra                              
param Pgb{BAR};                 #Potencia Activa generada en la barra
param Qgb{BAR};                 #Potencia Reactiva generada en la barra                               
param gshb{BAR};                #Conductancia shunt en p.u. conectada a la barra
param bshb{BAR};                #Susceptancia shunt en p.u. conectada a la barra
param wshm{BAR};                #Estado del equipo de compensacion reactiva 1 si esta conectado y 0 si no esta conectado                               
param Pgm{BAR};                 #Potencia Activa minima generada en la barra
param PgM{BAR};                 #Potencia Activa maxima generada en la barra
param Qgm{BAR};                 #Potencia Reactiva minima generada en la barra
param QgM{BAR};                 #Potencia Reactiva maxima generada en la barra
param c2{BAR};                  #coeficiente de la variable cuadratica de la funcion costo del generador
param c1{BAR};                  #coeficiente de la variable lineal de la funcion costo del generador
param c0{BAR};                  #coeficiente independiente del costo de generacion

# Parametros de los ramos
param r{RAM};					#Resistencia de la linea en p.u.
param x{RAM};                   #Reactancia de la linea en p.u.
param bsh{RAM};                 #Susceptancia shunt de la linea en p.u.
param fi{RAM};                  #Angulo de desfase de transformadores cambiadores de fase en grados
param angmin{RAM};              #Angulo minimo permitido de desfase angular de la linea
param angmax{RAM};              #Angulo maximo permitido de desfase angular de la linea
param a{RAM};                   #Relacion de transformacion del transformador
param SMAX{RAM};                #Capacidad maxima de flujo de potencia aparente de la linea de transmision
param ntap{RAM};                #Numero de posiciones discretos del Tap en el transformador
param rMtap{RAM};               #Relacion maxima de transformacion del transformador
param rmtap{RAM};               #Relacion minima de transformacion del transformador
 
# parametros generales
param Sbase;					#Potencia base del sistema
param Vnom;                     #Tension nominal p.u. del sistema
param Vmin;                     #Tension nominal minimo en p.u. de operacion permitido para las barras
param Vmax;                     #Tension nominal maximo en p.u. de operacion permitido para las barras

# Parametros adicionales para el FPO_IC
param g{RAM};					#Conductancia de la linea de transmision en p.u.
param b{RAM};                   #Susceptancia de la linea de transmision en p.u.
param reg{RAM};                 #Regulacion del tap del transformador en la linea

#Declaracion de variables de ramos
var nt{RAM} integer; 			#Tap inicial del tranformador en la linea
var Ire{RAM};					#Parte Real de la corriente en la linea
var Iim{RAM};					#Parte Imaginaria de la corriente en la linea

# Declaracion de las variables en las barras
var w{BAR} binary; 				#Estado del Equipo de compensacion reactiva en la barra
var Pger{BAR} >= 0;				#Potencia activa generada en la barra
var Qger{BAR};                  #Potencia reactiva generada en la barra
var Vsqr{BAR} >= 0;				#Raiz cuadrada de la tension en la barra
var ISre{BAR};					#Parte real de la Corriente total de la barra
var ISim{BAR};					#Parte imaginaria de la Corriente total de la barra
var IDre{BAR};					#Parte real de la Corriente demandada de la barra
var IDim{BAR};					#Parte imaginaria de la Corriente demandada de la barra
var Vre{BAR};					#Parte real de la tensión en la barra
var Vim{BAR};                   #Parte imaginaria de la tensión en la barra

#################################################
# Definicion de la funcion objetivo

minimize custo_generacion:
	sum {i in BAR: tipo[i] <> 0}(c2[i]*(Pger[i]*Sbase)^2 + c1[i]*Pger[i]*Sbase + c0[i]);

##################################################
# Definicion de las restriciones

#Balance de corriente
subject to cre_bal {i in BAR}:
	sum {(k,j,i) in RAM}((1/(1 + reg[k,j,i] * (nt[k,j,i]/ntap[k,j,i])))*Ire[k,j,i] + bsh[k,j,i]*Vim[i]) - sum{(k,i,j) in RAM}(Ire[k,i,j] - bsh[k,i,j]*Vim[i]) + ISre[i] + bshb[i]*Vim[i]*w[i] = IDre[i];

subject to cim_bal {i in BAR}:
	sum {(k,j,i) in RAM}((1/(1 + reg[k,j,i] * (nt[k,j,i]/ntap[k,j,i])))*Iim[k,j,i] - bsh[k,j,i]*Vre[i]) - sum{(k,i,j) in RAM}(Iim[k,i,j] + bsh[k,i,j]*Vre[i]) + ISim[i] - bshb[i]*Vre[i]*w[i] = IDim[i];

#Caida de tension
subject to vre_dr {(k,i,j) in RAM}:
	(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*Vre[i] - Vre[j] = (1/(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])))*(r[k,i,j]*Ire[k,i,j] - x[k,i,j]*Iim[k,i,j]);

subject to vim_dr {(k,i,j) in RAM}:
	(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*Vim[i] - Vim[j] = (1/(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])))*(x[k,i,j]*Ire[k,i,j] + r[k,i,j]*Iim[k,i,j]);

#Corriente de carga
subject to cre_ldnl {i in BAR: Pd[i] <> 0 || Qd[i] <> 0}:
	IDre[i] = (Pd[i]*Vre[i] + Qd[i]*Vim[i])/Vsqr[i];

subject to cim_ldnl {i in BAR: Pd[i] <> 0 || Qd[i] <> 0}:
	IDim[i] = (Pd[i]*Vim[i] - Qd[i]*Vre[i])/Vsqr[i];

subject to calculo_potencia_activa_generada {i in BAR}:
	Pger[i] = Vre[i]*ISre[i] + Vim[i]*ISim[i];

subject to calculo_potencia_reactiva_generada {i in BAR}:
	Qger[i] = Vim[i]*ISre[i] - Vre[i]*ISim[i];

subject to calculo_magnitude_tension {i in BAR}:
	Vsqr[i] = Vre[i]^2 + Vim[i]^2;

######################################################
# Definicion de los limites

#Limites de generacion
subject to limite_potencia_activa_generada {i in BAR}:
	Pgm[i] <= Pger[i] <= PgM[i];

subject to limite_potencia_reactiva_generada {i in BAR}:
	Qgm[i] <= Qger[i] <= QgM[i];

#Limites de tension
subject to vl_lim {i in BAR}:
	Vmin^2 <= Vsqr[i] <= Vmax^2;

#Limites de las lineas
subject to limite_flujo_de {(k,i,j) in RAM}:
	Vsqr[i]*((Ire[k,i,j] - bsh[k,i,j]*Vim[i])^2 + (Iim[k,i,j] + bsh[k,i,j]*Vre[i])^2) <= SMAX[k,i,j]^2;

subject to limite_flujo_para {(k,i,j) in RAM}:
	Vsqr[j]*(((1/(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])))*Ire[k,i,j] + bsh[k,i,j]*Vim[i])^2 + ((1/(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])))*Iim[k,i,j] - bsh[k,i,j]*Vre[i])^2) <= SMAX[k,i,j]^2;

# limite de nt
subject to limite_nt {(k,i,j) in RAM}:
	- ntap[k,i,j] <= nt[k,i,j] <= ntap[k,i,j];

# limite de w
subject to limite_w {i in BAR}:
	0 <= w[i] <= wshm[i];