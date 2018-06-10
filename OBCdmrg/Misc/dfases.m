
function dfases(filehm,oya,S1,S0,celdaU,Nrp)
% DFASES - Construye el diagrama de fase cuantico

% Esto carga la variable hm (magnetizacion vs campo h)
load(filehm)

% Encontrar los valores de magnetizacion persistente
lD = size(hm,3); pl = cell(1,lD);
for k = 1:lD
    x = [];
    y = hm(:,2,k);
    while (~isempty(y))
        mtch = find(abs(y(1)-y)<1e-12);
        if numel(mtch) > 1
            x(end+1) = y(1);
        end
        y(mtch) = [];
    end
    pl{k} = x;
end

% Seleccionar una lista unica x de estos valores
x = [];
y = [pl{:}]; clear pl
while (~isempty(y))
  miny = min(y); x(end+1) = miny;
  y = y(abs(miny-y)>1e-12);
end

% Eliminar errores de precision -> Comparar con valores permitidos
vecP  = repmat(celdaU,[1 Nrp]);
SZmax = S0*numel(find(vecP==0)) + S1*numel(find(vecP==1));
Mag   = (-SZmax:1:SZmax)./(2*Nrp); ind = []; 
for n = 1:length(x)
    indn = find(abs(x(n)-Mag)<1e-12);
    if ~isempty(indn), ind(end+1) = indn; end
end
m = Mag(ind); lm = length(m); % Comparar x y m

% Si se quieren filtrar plateaus de acuerdo al criterio OYA
if oya, [m,lm] = criterioOYA(m,S1,S0); end

tol = oya*0.1 + (1-oya)*1e-12;
% Eliminar errores
for k = 1:lD
    for n = 1:lm
        ind = find(abs(hm(:,2,k)-m(n))<tol);
        if isempty(ind), continue, end
        hm(ind,2,k) = m(n);
    end
end

% Campos criticos
hc = cell(1,lm-1);
for k = 1:lD
    for n = 1:lm-1
        ind = find(hm(:,2,k)==m(n+1)); 
        if isempty(ind), continue, end
        hc{n} = [hc{n};[0.1*(k-1),hm(ind(1),1,k)]];
    end
end

% Simbolos para dibujar
%sb = ['+','x','*','^','.','d','s','v','<','>','o','p','h'];
%sb = ['v','^','s','.','d','s','v','<','>','o','p','h'];
sb = ['o','^','s','.','d','s','v','<','>','o','p','h'];
% Color del simbolo
cs = ['b','g','r','k','m','c','y','w'];
% Unir con linea de tipo
lt = ['-',':','-.','--','(none)'];

% Dibujar diagrama de fases
hold on
for n = 1:lm-1
    %lin = strcat(lt(1),cs(4),sb(n));
    lin = strcat(cs(4),sb(n));
    plot(hc{n}(:,1),hc{n}(:,2),lin,'MarkerFaceColor',cs(4))
end
   
% % Utilizar Latex para simbolos matematicos en texto
% Interpreter = {'Interpreter','latex'};
% % Tamaño de la fuente en texto
% FontSize = {'FontSize',12};
% % Titulo de la figura
% Title = ['Diagrama de fase cuántico, Ising (', ...
%     num2str(S1),',',num2str(S0),')'];
% title(Title,FontSize{:})
% xlabel('D',FontSize{:}), ylabel('h',FontSize{:})

%=========================================================================
function [m,lm] = criterioOYA(x,S1,S0)
% CRITERIOOYA - Aplica el criterio OYA para ferrimagneticos

% Aplicar condicion de OYA: (S0+S1)*(1-m/msat) = entero
% Los plateaus que no cumplan esto se eliminan
msat = max(x);
S0mS1 = S0+S1;
ind = [];
for n = 1:length(x);
    Noya = S0mS1*(1-x(n)/msat);
    if Noya == fix(Noya)
        ind(end+1) = n;
    end
end
m = x(ind); lm = length(m);
%=========================================================================
                            