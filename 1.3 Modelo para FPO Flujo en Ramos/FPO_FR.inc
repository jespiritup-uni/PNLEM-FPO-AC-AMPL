# Codigo para imprimir los modelos
printf "\n=====================================================";
printf "\nRESUMEN DEL MODELO FLUJO DE POTENCIA OPTIMO BASADO DE FLUJO EN RAMOS";
printf "\n=====================================================";

printf "\nCosto de Generacion (US$)         : ";
printf "%10.6f",custo_generacion;
printf "\n";
printf "\nPerdidas            (MW)         :";
printf "%10.4f",sum {(k,i,j) in RAM}(Sbase*(r[k,i,j]*Iqdr[k,i,j]));
printf "\n";

printf "\nTENSION Y ANGULO EN LA BARRA \n";
printf "barra     V[pu]    theta[grau]   Pger[MW]     Qger[MW]      w                Balance activo                          Balance reactivo\n";
for {i in BAR} {
	printf "%5d",i;
	printf "%12.4f %12.4f %12.4f %12.4f %12.4f %25.20f %25.20f",
		sqrt(Vqdr[i]),180*theta[i]/3.14159265359,
		Sbase*Pger[i],Sbase*Qger[i],
		w[i],
		Pger[i] - Pd[i] + sum {(k,j,i) in RAM}(P[k,j,i]) - sum {(k,i,j) in RAM}(P[k,i,j] + r[k,i,j]*Iqdr[k,i,j]) - Vqdr[i] * gshb[i],
		Qger[i] - Qd[i] + sum {(k,j,i) in RAM}(Q[k,j,i] + bsh[k,j,i]*Vqdr[i]) - sum {(k,i,j) in RAM}(Q[k,i,j] - bsh[k,i,j]*Vqdr[i] + x[k,i,j]*Iqdr[k,i,j]) + Vqdr[i]*bshb[i]*w[i];
	printf "\n";
}
printf "\n";

printf "\nFLUJO DE POTENCIA ACTIVA Y REACTIVA, PERDIDAS\n";
printf "     i      j        Pij[MW]        Pji[MW]       Qij[MVAr]      Qji[MVAr]     Pls[MW]     Qls[MVAr]       Tap       Pos_tap       Ang[graus]\n";
for {(k,i,j) in RAM} {
	printf "%5d %5d",i,j;
	printf "%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f",
	Sbase*(P[k,i,j]+r[k,i,j]*Iqdr[k,i,j]),-Sbase*(P[k,i,j]),
	Sbase*(Q[k,i,j]-bsh[k,i,j]*Vqdr[i]+x[k,i,j]*Iqdr[k,i,j]),-Sbase*(Q[k,i,j]+bsh[k,i,j]*Vqdr[j]),
	Sbase*(r[k,i,j]*Iqdr[k,i,j]),Sbase*(x[k,i,j]*Iqdr[k,i,j] -bsh[k,i,j]*Vqdr[i] - bsh[k,i,j]*Vqdr[j]),
	(1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])),nt[k,i,j],
	180*theta[i]/3.14159265359-180*theta[j]/3.14159265359;
	printf "\n";
}
printf "\n";

printf "\nTotal CPU time:%10.2fs", _ampl_time + _total_solve_time + _total_shell_time;
printf "\nElapsed time :%10.2fs\n", _solve_elapsed_time;
printf "=========================================================================================================================================================================\n";