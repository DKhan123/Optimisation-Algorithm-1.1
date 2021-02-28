% Determine Call sign ID needed for thermodynamic parameters 
function GEMId=CallSign(Thermomodel)
ThermoRead=fileread(Thermomodel);
ThermoIdn=extractBetween(ThermoRead,'order):','See');
ThermoIdn=split(strip(extractBetween(ThermoIdn,'|','|')));
ThermoIdn=strrep(ThermoIdn,"'",'prime');
ThermoIdn=strrep(ThermoIdn,"-",'');
ThermoRead=strjoin(splitlines(strjoin(extractBetween(ThermoRead,')','|'))));
ThermoRead=string(strip(split(ThermoRead,'=')));
ThermoRead=unique(strip(eraseBetween(ThermoRead,1,' ')))';
ThermoRead(cellfun(@length, ThermoRead)~=2)=[];
ThermoRead(1)= '';
PerpId=["f0" "n" "V0"  "k0" "k0prime" "td"  "gam0"  "q"  "etaS0" "Sconf" "g0" "g0prime"];
GEMId=[PerpId;ThermoRead(1:length(ThermoIdn))+' = '];

end

