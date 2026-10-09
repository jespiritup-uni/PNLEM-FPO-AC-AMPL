# Declaracion de conjuntos
set BAR;										# Conjunto de barras
set RAM within 1..500 cross BAR cross BAR;		# Conjunto de ramos

# Parametros de las barras
param tipo{BAR};			#Tipo de barra (1=PQ o carga, 2=PV o generacion, 3=Slack o referencia)
param Vb{BAR};              #Magnitud de tensión inicial o de operacion en la barra en p.u.
param angi{BAR};            #Angulo de la tension inical en la barra en grados
param Pd{BAR};              #Potencia Activa demandada en la barra
param Qd{BAR};              #Potencia Reactiva demandada en la barra                           
param Pgb{BAR};             #Potencia Activa generada en la barra
param Qgb{BAR};             #Potencia Reactiva generada en la barra                           
param gshb{BAR};            #Conductancia shunt en p.u. conectada a la barra
param bshb{BAR};            #Susceptancia shunt en p.u. conectada a la barra
param wshm{BAR};            #Estado del equipo de compensacion reactiva 1 si esta conectado y 0 si no esta conectado                          
param Pgm{BAR};             #Potencia Activa minima generada en la barra
param PgM{BAR};             #Potencia Activa maxima generada en la barra
param Qgm{BAR};             #Potencia Reactiva minima generada en la barra
param QgM{BAR};             #Potencia Reactiva maxima generada en la barra
param c2{BAR};              #coeficiente de la variable cuadratica de la funcion costo del generador
param c1{BAR};              #coeficiente de la variable lineal de la funcion costo del generador
param c0{BAR};              #coeficiente independiente del costo de generacion

# Parametros de los ramos
param r{RAM};               #Resistencia de la linea en p.u.
param x{RAM};               #Reactancia de la linea en p.u.
param bsh{RAM};             #Susceptancia shunt de la linea en p.u.
param fi{RAM};              #Angulo de desfase de transformadores cambiadores de fase en grados
param angmin{RAM};          #Angulo minimo permitido de desfase angular de la linea
param angmax{RAM};          #Angulo maximo permitido de desfase angular de la linea
param a{RAM};               #Relacion de transformacion del transformador
param SMAX{RAM};            #Capacidad maxima de flujo de potencia aparente de la linea de transmision
param ntap{RAM};            #Numero de posiciones discretos del Tap en el transformador
param rMtap{RAM};           #Relacion maxima de transformacion del transformador
param rmtap{RAM};           #Relacion minima de trasnformacion del transformador

# parametros generales
param Sbase;				#Potencia base del sistema
param Vnom;                 #Tension nominal p.u. del sistema
param Vmin;                 #Tension nominal minimo en p.u. de operacion permitido para las barras
param Vmax;                 #Tension nominal maximo en p.u. de operacion permitido para las barras

# Parametros adicionales para el FPO_CP
param g{RAM};				#Conductancia de la linea de transmision en p.u.
param b{RAM};               #Susceptancia de la linea de transmision en p.u.
param reg{RAM};             #Regulacion del tap del transformador en la linea

#Declaracion de variables de ramos
var nt{RAM} integer;        #Tap inicial del tranformador en la linea
var Ppara{RAM};				#Potencia activa de salida de la linea de transmision
var Qpara{RAM};             #Potencia reactiva de salida de la linea de transmision
var Pde{RAM};               #Potencia activa de entrada de la linea de transmision
var Qde{RAM};               #Potencia reactiva de entrada de la linea de transmision

# Declaracion de las variables en las barras
var VR{BAR};   				#Parte real de la tensión en la barra
var VT{BAR};   				#Parte imaginaria de la tensión en la barra

	# Variables Auxiliares: tensión del lado emisor después de tap a y desfase fi
var Vsq{BAR};				#Raiz cuadrada de la tension en la barra 
var VRp{RAM};         	 	#Parte real de la tension despues del transformador en la linea, VR' = (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*(VR*cos(fi) - VT*sin(fi))
var VTp{RAM};          		#Parte imaginaria de la tension despues del transformador en la linea, VT' = (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*(VR*sin(fi) + VT*cos(fi))
var Pger{BAR};				#Potencia activa generada en la barra
var Qger{BAR};              #Potencia reactiva generada en la barra
var gP{BAR} >= 0;			#Potencia activa en la barra
var gQ{BAR} >= 0;           #Potencia reactiva en la barra
var w{BAR} binary;			#Estado del equipo de compensacion reactiva en la barra

#################################################
# Definicion de la funcion objetivo

minimize custo_generacion:
	sum {i in BAR: tipo[i] <> 0}(c2[i]*(Pger[i]*Sbase)^2 + c1[i]*Pger[i]*Sbase + c0[i]);

##################################################
# Definicion auxiliares: VRp y VTp

subject to def_VRp {(k,i,j) in RAM}:
    VRp[k,i,j] = (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*( VR[i]*cos(fi[k,i,j]) - VT[i]*sin(fi[k,i,j]) );

subject to def_VTp {(k,i,j) in RAM}:
    VTp[k,i,j] = (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*( VR[i]*sin(fi[k,i,j]) + VT[i]*cos(fi[k,i,j]) );

##################################################
# Definicion de las restriciones

# Balance de potencia
subject to balance_P {i in BAR}:
    Pger[i] - Pd[i]
    - sum {(k,j,i) in RAM} ( Ppara[k,j,i] )
    - sum {(k,i,j) in RAM} ( Pde[k,i,j] )
    - gshb[i]*(Vsq[i]) = 0;

subject to balance_Q {i in BAR}:
    Qger[i] - Qd[i]
    - sum {(k,j,i) in RAM} ( Qpara[k,j,i] )
    - sum {(k,i,j) in RAM} ( Qde[k,i,j] )
    + bshb[i]*(Vsq[i])*w[i] = 0;

# Potencia activa desde i hacia j (Pde)
subject to flujo_P_de {(k,i,j) in RAM}:
    Pde[k,i,j] =
        g[k,i,j]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*(Vsq[i])
        - (
            g[k,i,j]*( VRp[k,i,j]*VR[j] + VTp[k,i,j]*VT[j] )
            + b[k,i,j]*( VTp[k,i,j]*VR[j] - VRp[k,i,j]*VT[j] )
          );

# Potencia activa desde j hacia i (Ppara)
subject to flujo_P_para {(k,i,j) in RAM}:
    Ppara[k,i,j] =
        g[k,i,j]*(VR[j]^2 + VT[j]^2)
        - (
            g[k,i,j]*( VRp[k,i,j]*VR[j] + VTp[k,i,j]*VT[j] )
            - b[k,i,j]*( VTp[k,i,j]*VR[j] - VRp[k,i,j]*VT[j] )
          );

# Potencia reactiva desde i hacia j (Qde)
subject to flujo_Q_de {(k,i,j) in RAM}:
    Qde[k,i,j] =
        -(b[k,i,j] + bsh[k,i,j])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*(Vsq[i])
        + (
            b[k,i,j]*( VRp[k,i,j]*VR[j] + VTp[k,i,j]*VT[j] )
            - g[k,i,j]*( VTp[k,i,j]*VR[j] - VRp[k,i,j]*VT[j] )
          );

# Potencia reactiva desde j hacia i (Qpara)
subject to flujo_Q_para {(k,i,j) in RAM}:
    Qpara[k,i,j] =
        -(b[k,i,j] + bsh[k,i,j])*(VR[j]^2 + VT[j]^2)
        + (
            b[k,i,j]*( VRp[k,i,j]*VR[j] + VTp[k,i,j]*VT[j] )
            + g[k,i,j]*( VTp[k,i,j]*VR[j] - VRp[k,i,j]*VT[j] )
          );
		  
######################################################
# Definicion de los limites

# limite nt
subject to limite_nt {(k,i,j) in RAM}:
	- ntap[k,i,j] <= nt[k,i,j] <= ntap[k,i,j];

# limite do w
subject to limite_w {i in BAR}:
	0 <= w[i] <= wshm[i];

# limite do tension
subject to limite_tension_1 {i in BAR}:
	Vsq[i] = VR[i]^2 + VT[i]^2;

subject to limite_tension {i in BAR}:
    Vmin^2 <= Vsq[i] <= Vmax^2;

# Limite do potencia P e Q
subject to limite_potencia_generada_activa {i in BAR}:
	Pgm[i] <= Pger[i] <= PgM[i];

subject to limite_potencia_generada_reactiva {i in BAR}:
	Qgm[i] <= Qger[i] <= QgM[i];

# cuadrado de flujo de potencia P e Q
subject to cuadrado_flujo_potencia_de {(k,i,j) in RAM}:
    Pde[k,i,j]^2 + Qde[k,i,j]^2 <= SMAX[k,i,j]^2;

subject to cuadrado_flujo_potencia_para {(k,i,j) in RAM}:
    Ppara[k,i,j]^2 + Qpara[k,i,j]^2 <= SMAX[k,i,j]^2;
