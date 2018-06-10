
load hm

% Encontrar los valores de magnetizacion persistente
lD = size(hm,3); pl = cell(1,lD);
for k = 1:lD
    x = [];
    y = hm(:,2,k);
    while (~isempty(y))
        mtch = find(y(1) == y);
        if numel(mtch) > 1
            x(end+1) = y(1);
        end
        y(mtch) = [];
    end
    x(abs(x)<eps)=0;
    pl{k} = x;
end

% Seleccionar una lista unica de estos valores
x = [];
y = [pl{:}];
while (~isempty(y))
  x(end+1) = min(y);
  y = y(y~=min(y));
end

% Aplicar condicion de OYA: (S0+S1)*(1-m/msat) = entero
% Los plateaus que no cumplan esto se eliminan
S0 = 1.5; S1 = 1.5;
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

% Eliminar plateaus falsos
hmag = hm;
for k = 1:lD
    for j = 1:size(hmag,1)
        for n = 1:lm 
            if abs(hmag(j,2,k)-m(n)) < 0.1
                hmag(j,2,k) = m(n);
            else
                continue
            end
        end
    end
end

% Campos criticos
hc = cell(1,lm-1);
for k = 1:lD
    for n = 1:lm-1
        ind = find(hmag(:,2,k)== m(n+1)); 
        if isempty(ind), continue, end
        hc{n} = [hc{n};[0.1*(k-1),hmag(ind(1),1,k)]];
    end
end

% Simbolos para dibujar
sb = ['o','>','v','^'];
% Color del simbolo
cs = 'k';
% Unir con linea de tipo
lt = '-';

% Dibujar diagrama de fases
hold on
for n = 1:lm-1
    lin = strcat(lt,cs,sb(n));
    plot(hc{n}(:,1),hc{n}(:,2),lin,'MarkerFaceColor',cs,'MarkerSize',5)
end
