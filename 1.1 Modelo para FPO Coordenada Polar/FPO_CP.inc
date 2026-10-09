# Codigo para imprimir los modelos
printf "\n==========================================================";
printf "\nRESUMEN DEL MODELO FLUJO DE POTENCIA OPTIMA BASADO EN COORDENADAS POLARES";
printf "\n==========================================================";

printf "\nCosto de Generacion (US$/h)         : ";
printf "%10.6f",custo_generacion;
printf "\n";
printf "\nPerdidas            (MW)         :";
printf "%10.4f",sum {(k,i,j) in RAM}(Sbase*(Pde[k,i,j]+Ppara[k,i,j]));
printf "\n";

printf "\nTENSION Y ANGULO EN LA BARRA \n";
printf "barra     V[pu]    theta[grau]   Pger[MW]     Qger[MW]      w                Balance activo                          Balance reactivo\n";
for {i in BAR} {
	printf "%5d",i;
	printf "%12.4f %12.4f %12.4f %12.4f %12.4f %25.20f %25.20f",
		V[i],180*th[i]/3.14159265359,
		Sbase*Pger[i],Sbase*Qger[i],
		w[i],
		Pger[i] - Pd[i] - sum {(k,j,i) in RAM}(Ppara[k,j,i]) - sum {(k,i,j) in RAM}(Pde[k,i,j]) - V[i]^2*gshb[i],
		Qger[i] - Qd[i] - sum {(k,j,i) in RAM}(Qpara[k,j,i]) - sum {(k,i,j) in RAM}(Qde[k,i,j]) + V[i]^2*bshb[i]*w[i];
	printf "\n";
}
printf "\n";

printf "\nFLUJO DE POTENCIA ACTIVA Y REACTIVA, PERDIDAS E ANGULOS\n";
printf "     i      j        Pij[MW]        Pji[MW]       Qij[MVAr]      Qji[MVAr]     Pls[MW]     Qls[MVAr]       Tap       Pos_tap       Ang[graus]\n";
for {(k,i,j) in RAM} {
	printf "%5d %5d",i,j;
	printf "%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %30.20f %30.20f %30.20f %30.20f",
	Sbase*Pde[k,i,j],Sbase*Ppara[k,i,j],
	Sbase*Qde[k,i,j],Sbase*Qpara[k,i,j],
	Sbase*(Pde[k,i,j]+Ppara[k,i,j]),Sbase*(Qde[k,i,j]+Qpara[k,i,j]),
	(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])),nt[k,i,j],
	(th[i] - th[j]) * 180 / 3.14159265359, 
	Pde[k,i,j] - (g[k,i,j]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]))),
	Ppara[k,i,j] - (g[k,i,j]*V[j]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]))),
	Qde[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j]))),
	Qpara[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*V[j]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])))
	;
	printf "\n";
}
printf "\n";

printf "\nBalance activo\n";
for {i in BAR} {
	if (abs(Pger[i] - Pd[i] - sum {(k,j,i) in RAM}(Ppara[k,j,i]) - sum {(k,i,j) in RAM}(Pde[k,i,j]) - V[i]^2*gshb[i]) > 1e-9) then{
	printf "%5d %25.20f\n",i, Pger[i] - Pd[i] - sum {(k,j,i) in RAM}(Ppara[k,j,i]) - sum {(k,i,j) in RAM}(Pde[k,i,j]) - V[i]^2*gshb[i];
	}
}

printf "\nBalance reactivo\n";
for {i in BAR} {
	if (abs(Qger[i] - Qd[i] - sum {(k,j,i) in RAM}(Qpara[k,j,i]) - sum {(k,i,j) in RAM}(Qde[k,i,j]) + V[i]^2*bshb[i]*w[i]) > 1e-9) then{
	printf "%5d %25.20f\n",i, Qger[i] - Qd[i] - sum {(k,j,i) in RAM}(Qpara[k,j,i]) - sum {(k,i,j) in RAM}(Qde[k,i,j]) + V[i]^2*bshb[i]*w[i];
	}
}

printf "\nPde\n";
for {(k,i,j) in RAM} {
	if (abs(Pde[k,i,j] - (g[k,i,j]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])))) > 1e-7) then{
	printf "%5d %5d %25.20f\n",i,j, Pde[k,i,j] - (g[k,i,j]*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])));
	}
}

printf "\nPpara\n";
for {(k,i,j) in RAM} {
	if (abs(Ppara[k,i,j] - (g[k,i,j]*V[j]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])))) > 1e-7) then{
	printf "%5d %5d %25.20f\n",i,j, Ppara[k,i,j] - (g[k,i,j]*V[j]^2 - (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(g[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - b[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])));
	}
}

printf "\nQde\n";
for {(k,i,j) in RAM} {
	if (abs(Qde[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])))) > 1e-7) then{
	printf "%5d %5d %25.20f\n",i,j, Qde[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))^2*V[i]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) - g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])));
	}
}

printf "\nQpara\n";
for {(k,i,j) in RAM} {
	if (abs(Qpara[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*V[j]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])))) > 1e-7) then{
	printf "%5d %5d %25.20f\n",i,j, Qpara[k,i,j] - (-(b[k,i,j] + bsh[k,i,j])*V[j]^2 + (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j]))*V[i]*V[j]*(b[k,i,j]*cos(th[i] - th[j] + fi[k,i,j]) + g[k,i,j]*sin(th[i] - th[j] + fi[k,i,j])));
	}
}

printf "\nViolacion de tension\n";
for {i in BAR} {
	if (V[i] > Vmax) then{
	printf "%5d %25.20f\n",i, V[i]-Vmax;
	}

	if (V[i] < Vmin) then{
	printf "%5d %25.20f\n",i, Vmin-V[i];
	}
}

printf "\nViolacion de generacion\n";
for {i in BAR} {
	if (Pger[i] > PgM[i]) then{
	printf "%5d %25.20f\n",i, Sbase*(Pger[i]-PgM[i]);
	}

	if (Qger[i] > QgM[i]) then{
	printf "%5d %25.20f\n",i, Sbase*(Qger[i]-QgM[i]);
	}

	if (Pger[i] < Pgm[i]) then{
	printf "%5d %25.20f\n",i, Sbase*(Pgm[i]-Pger[i]);
	}

	if (Qger[i] < Qgm[i]) then{
	printf "%5d %25.20f\n",i, Sbase*(Qgm[i]-Qger[i]);
	}
}

printf "\nViolacion de capacidad de la linea\n";
for {(k,i,j) in RAM} {
	if (Pde[k,i,j]^2 + Qde[k,i,j]^2 > SMAX[k,i,j]^2) then{
	printf "%5d %5d %25.20f\n",i,j, ((Sbase*Pde[k,i,j])^2 + (Sbase*Qde[k,i,j])^2) - (Sbase*SMAX[k,i,j])^2;
	}

	if (Ppara[k,i,j]^2 + Qpara[k,i,j]^2 > SMAX[k,i,j]^2) then{
	printf "%5d %5d %25.20f\n",i,j, ((Sbase*Ppara[k,i,j])^2 + (Sbase*Qpara[k,i,j])^2) - (Sbase*SMAX[k,i,j])^2;
	}
}

printf "\nTotal CPU time:%10.2fs", _ampl_time + _total_solve_time + _total_shell_time;
printf "\nElapsed time :%10.2fs\n", _solve_elapsed_time;
printf "=========================================================================================================================================================================\n";
