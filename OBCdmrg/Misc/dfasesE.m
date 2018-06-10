
clear 
load he

lD = size(he,3); lh = size(he,1);

M = zeros(lh,lD);
for k = 1:lD
    for n = 1:lh-1
        M(n,k) = (he(n,2,k)-he(n+1,2,k))/(he(n+1,1,k)-he(n,1,k));
    end
end
M(lh,:) = M(lh-1,:);

% Encontrar los valores de magnetizacion persistente
pl = cell(1,lD);
for k = 1:lD
    x = [];
    y = M(:,k);
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
S1 = 1; S0 = 1.5; celdaU = [1 0]; Nrp = 30;
vecP  = repmat(celdaU,[1 Nrp]);
SZmax = S0*numel(find(vecP==0)) + S1*numel(find(vecP==1));
Mag   = (-SZmax:1:SZmax)./(2*Nrp); ind = []; 
for n = 1:length(x)
    indn = find(abs(x(n)-Mag)<1e-12);
    if ~isempty(indn), ind(end+1) = indn; end
end
m = Mag(ind); lm = length(m); % Comparar x y m

tol = 1e-12;
% Eliminar errores
for k = 1:lD
    for n = 1:lm
        ind = find(abs(M(:,k)-m(n))<tol);
        if isempty(ind), continue, end
        M(ind,k) = m(n);
    end
end

% Campos criticos
hc = cell(1,lm-1);
for k = 1:lD
    for n = 1:lm-1
        ind = find(M(:,k)==m(n+1)); 
        if isempty(ind), continue, end
        hc{n} = [hc{n};[0.1*(k-1),he(ind(1),1,k)]];
    end
end

% Simbolos para dibujar
sb = ['o','d','v','^','<','>','s','.','p','h','+','*','x'];
% Color del simbolo
cs = ['b','g','r','k','m','c','y','w'];
% Unir con linea de tipo
lt = ['-',':','-.','--','(none)'];


% Dibujar diagrama de fases
hold on
for n = 1:lm-1
    lin = strcat(lt(1),cs(n),sb(1));
    plot(hc{n}(:,1),hc{n}(:,2),lin,'MarkerFaceColor',cs(n))
end