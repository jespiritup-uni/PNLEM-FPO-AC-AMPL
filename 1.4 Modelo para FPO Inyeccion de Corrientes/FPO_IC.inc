# Codigo para imprimir los modelos
printf "\n================================================================";
printf "\nRESUMEN DEL MODELO FLUJO DE POTENCIA OPTIMA BASADO EN INYECCION DE CORRIENTES";
printf "\n==========================================================";

printf "\nCosto de Generacion (US$/h)         : ";
printf "%10.6f\n",custo_generacion;
printf "\n";
printf "\nPerdidas            (MW)         :";
printf "%10.4f",sum {(k,i,j) in RAM}(Sbase*r[k,i,j]*(Ire[k,i,j]^2 + Iim[k,i,j]^2));
printf "\n";

printf "TENSION Y ANGULO EN LA BARRA\n";
printf "barra     V[pu]    theta[grau]   Pger[MW]     Qger[MW]      w                Balance activo                          Balance reactivo\n";
for {i in BAR} {
    printf "%5d",i;
	printf "%12.4f %12.4f %12.4f %12.4f %12.4f %25.20f %25.20f\n",
    sqrt(Vre[i]^2 + Vim[i]^2),
    atan2(Vim[i],Vre[i])*180/(4*atan(1)),
    Pger[i] * Sbase,
    Qger[i] * Sbase,
	w[i],
    cre_bal[i].slack,
    cim_bal[i].slack;
}
printf "\n";

printf "FLUJO DE POTENCIA ACTIVA Y REACTIVA, PERDIDAS Y ANGULOS\n";
printf "     i      j        Pij[MW]        Pji[MW]       Qij[MVAr]      Qji[MVAr]     Pls[MW]     Qls[MVAr]       Tap       Pos_tap       Ang[graus]\n";
for {(k,i,j) in RAM} {
    printf "%5d %5d",i,j;
	printf "%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f\n",
    (Vre[i]*Ire[k,i,j] + Vim[i]*Iim[k,i,j]) * Sbase,
    (Vre[j]*(-Ire[k,i,j]) + Vim[j]*(-Iim[k,i,j])) * Sbase,
    (Vim[i]*Ire[k,i,j] - Vre[i]*Iim[k,i,j]) * Sbase,
    (Vim[j]*(-Ire[k,i,j]) - Vre[j]*(-Iim[k,i,j])) * Sbase,
    (Vre[i]*Ire[k,i,j] + Vim[i]*Iim[k,i,j] + Vre[j]*(-Ire[k,i,j]) + Vim[j]*(-Iim[k,i,j]) ) * Sbase,
    (Vim[i]*Ire[k,i,j] - Vre[i]*Iim[k,i,j] + Vim[j]*(-Ire[k,i,j]) - Vre[j]*(-Iim[k,i,j])) * Sbase,
    (1 + reg[k,i,j] * (nt[k,i,j]/ntap[k,i,j])),nt[k,i,j],
    (180/3.14159265359) * atan2(Vim[i], Vre[i]) - (180/3.14159265359) * atan2(Vim[j], Vre[j]);
}

printf "\n";
printf "\nTotal CPU time:%10.2fs", _ampl_time + _total_solve_time + _total_shell_time;
printf "\nElapsed time :%10.2fs\n", _solve_elapsed_time;
printf "=========================================================================================================================================================================\n";