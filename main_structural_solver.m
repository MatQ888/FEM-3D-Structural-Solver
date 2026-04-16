% Inicio del código, parámeros de la estructura
L = 1;
E = 100e9;
nu = 0.3;
a = L/20;
b = a/2; 
t = L/200;

P = 1000;
q0z = 500; %N/m
uload = L/100;
G = E/(2*(1+nu));


% Área (2 alas + 1 alma)
A = 2*(b*t) + (a-2*t)*t;

% Iz: Momento de inercia respecto al eje ze
% Es el rectángulo total (b*a^3/12) menos los dos huecos laterales
Iz = (b*a^3)/12 - ((b-t)*(a-2*t)^3)/12;

% Iy: Momento de inercia respecto al eje ye (EJE DÉBIL)
% Es la suma de las inercias de las alas y el alma respecto al centro
Iy = 2*(t*b^3)/12 + ((a-2*t)*t^3)/12;

% J: Módulo de torsión para perfil abierto
J = (1/3)*(2*b*t^3 + (a-2*t)*t^3);



%% Geometría
n1 = [L 0 0]; n2 = [L 0 (4/3)*L];
n3 = [L 0 (8/3)*L]; n4 = [L 0 4*L];
n5 = [L L 3.5*L]; n6 = [L 2*L 3*L];
n7 = [0 2*L 3*L];


% Direcciones de la base global en la matriz DIR
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]' ; DIR = [dirX, dirY, dirZ];
vaux1 = [0; -1; 0]; % Vector auxiliar de los elementos 1, 2 y 3
vaux2 = [0; 0; 1]; %vector auxiliar de los elementos 4, 5 y 6

%% Carga distribuida por secciones

qn2 = (1/3)*q0z;
qn3 = (2/3)*q0z;
qn4 = q0z;

% Elemento 1: n1 y n2
gdle1 = [1:12];
[Le1, Te1, vdire1] = Geom_viga3D(n1, n2, DIR, vaux1);
Ke_local1 = Keviga3D(E,G,Le1,A,J,Iy,Iz);
K_e1 = Te1'*Ke_local1*Te1;


%vector fuerzas nodales equivalente en ejes locales
% [u v w thx thy thz ...]
fneqe1 = -qn2*Le1*[0; 3/20; 0; 0; 0; Le1/30; 0; 7/20; 0; 0; 0; -Le1/20];
fneq_e1 = Te1'*fneqe1; % vector de fuerzas nodales equivalentes del elemento 1 en coordendadas globales



% Elemento 2: n2 y n3
gdle2 = [7:18];
[Le2, Te2, vdire2] = Geom_viga3D(n2, n3, DIR, vaux1);
Ke_local2 = Keviga3D(E,G,Le2,A,J,Iy,Iz);
K_e2 = Te2'*Ke_local2*Te2;

dq2 = qn3 - qn2;

% Parte Rectangular (Carga uniforme qn2)
% Posiciones Fy(2,8) y Mz(6,12)
f_rect = -qn2 * Le2 * [0; 1/2; 0; 0; 0; Le2/12; 0; 1/2; 0; 0; 0; -Le2/12];

% Parte Triangular (Carga extra dq2)
f_triang = -dq2 * Le2 * [0; 3/20; 0; 0; 0; Le2/30; 0; 7/20; 0; 0; 0; -Le2/20];

fneqe2 = f_rect + f_triang;
fneq_e2 = Te2' * fneqe2; % Siempre con traspuesta



% Elemento 3: n3 y n4
gdle3 = [13:24];
[Le3, Te3, vdire3] = Geom_viga3D(n3, n4, DIR, vaux1);
Ke_local3 = Keviga3D(E,G,Le3,A,J,Iy,Iz);
K_e3 = Te3'*Ke_local3*Te3;

f_rect3 = -(2*q0z/3) * Le3 * [0; 1/2; 0; 0; 0; Le3/12; 0; 1/2; 0; 0; 0; -Le3/12];
f_triang3 = -(q0z/3) * Le3 * [0; 3/20; 0; 0; 0; Le3/30; 0; 7/20; 0; 0; 0; -Le3/20];
fneq_e3 = Te3' * (f_rect3 + f_triang3);


% Elemento 4: n4 y n5
gdle4 = [19:30];
[Le4, Te4, vdire4] = Geom_viga3D(n4, n5, DIR, vaux2);
Ke_local4 = Keviga3D(E,G,Le4,A,J,Iy,Iz);
K_e4 = Te4'*Ke_local4*Te4;

fneq_e4 = zeros(12,1);

% Elemento 5: n5 y n6
gdle5 = [25:36];
[Le5, Te5, vdire5] = Geom_viga3D(n5, n6, DIR, vaux2);
Ke_local5 = Keviga3D(E,G,Le5,A,J,Iy,Iz);
K_e5 = Te5'*Ke_local5*Te5;

% Elemento 6: n6 y n7
gdle6 = [31:42];
[Le6, Te6, vdire6] = Geom_viga3D(n6, n7, DIR, vaux2);
Ke_local6 = Keviga3D(E,G,Le6,A,J,Iy,Iz);
K_e6 = Te6'*Ke_local6*Te6;


%% Ensamblaje
% Rigidez
K = zeros (42);
K (gdle1, gdle1) = K(gdle1, gdle1) + K_e1 % ensamblaje del ele 1
K (gdle2, gdle2) = K(gdle2, gdle2) + K_e2 % ensamblaje del ele 2
K (gdle3, gdle3) = K(gdle3, gdle3) + K_e3 % ensamblaje del ele 3
K (gdle4, gdle4) = K(gdle4, gdle4) + K_e4 % ensamblaje del ele 4
K (gdle5, gdle5) = K(gdle5, gdle5) + K_e5 % ensamblaje del ele 5
K (gdle6, gdle6) = K(gdle6, gdle6) + K_e6 % ensamblaje del ele 6

% Fuerzas
Fneq = zeros (42,1); % vector vacío de fuerzas nodales equivalentes del problema
Fneq (gdle1) = Fneq (gdle1) + fneq_e1; % ensamblaje del ele 1
Fneq (gdle2) = Fneq (gdle2) + fneq_e2; % ensamblaje del ele 2
Fneq (gdle3) = Fneq (gdle3) + fneq_e3; % ensamblaje del ele 3
Fneq (gdle4) = Fneq (gdle4) + fneq_e4; % ensamblaje del ele 4


%Condiciones de contorno
gdlT = [1:42]; % grados de libertad totales
gdlR = [1 2 3 37 38 39 40 41 42]; % grados de libertad restringidos
gdlU = [20]; %desplazamiento impuesto en el nodo 4 dirección 20
uload_val = L/100;
gdlL = setdiff(gdlT, [gdlR, gdlU]);

KLL = K(gdlL,gdlL);
KRR = K(gdlR, gdlR);
KLR = K (gdlL, gdlR);
KRL = KLR';


%% Desplazamientos
% Definimos uR (el vector de desplazamientos impuestos en las restricciones)
uR = zeros(length(gdlR), 1); 
% Si el desplazamiento uload está en la posición del gdlU dentro de gdlR,
% hay que asignarlo. Pero como tú definiste gdlU aparte:
% uR(donde_este_gdlU) = uload_val; 

% Vector de fuerzas totales (Equivalentes + Puntuales)
F_ext = Fneq; 
F_ext(27) = P; % Carga P en Nodo 5 (dirección Z global)

fL = F_ext(gdlL);
uL = KLL \ (fL - K(gdlL, gdlU) * uload_val);
Fneq_R = Fneq(gdlR);

% Reconstrucción del vector total U
U = zeros(42,1);
U(gdlL) = uL;
U(gdlU) = uload_val;
% U(gdlR) se queda a cero si son empotramientos

%% Calculo de reacciones
R_total = KRL*uL + KRR*uR - Fneq_R;
R1 = R_total (1:3);
R7 = R_total (4:9);


%% EXTRACCIÓN DE COMPONENTES
% Cada nodo i tiene: U(6*i-5)=Ux, U(6*i-4)=Uy, U(6*i-3)=Uz
Ux = U(1:6:42); % Desplazamientos en X de los 7 nodos
Uy = U(2:6:42); % Desplazamientos en Y de los 7 nodos
Uz = U(3:6:42); % Desplazamientos en Z de los 7 nodos

% Giros (opcional)
THx = U(4:6:42);
THy = U(5:6:42);
THz = U(6:6:42);

% Matriz de rigidez global (consultar página 37 de las diapositias)


%% Subrutinas

%% Viga 3D (beam element)
% Dirección de la base global en la matriz DIR
dirX = [1 0 0]'; dirY = [0 1 0]'; dirZ = [0 0 1]' ; DIR = [dirX, dirY, dirZ];
function [Le, Te, vdire] = Geom_viga3D(n1, n2, DIR, vaux)
    %[Le, Te] = Geom_viga3D([0 0 0], [1 0 0], DIR, [0;1;0])
    % n1,n2: coordenadas globales de los nodos locales 1 y 2
    % v_aux: vector auxiliar de orientación para el elemento - columna
    % base global de referencia
    dirX = DIR(:,1) ; dirY = DIR(:,2) ; dirZ = DIR(:,3);
    % deltas de longitud por coordenadas globales
    delta_x = n2(1)-n1(1) ; delta_y = n2(2)-n1(2); delta_z = n2(3)-n1(3);
    % Le longitud global del elemento
    Le = sqrt(delta_x^2 + delta_y^2 + delta_z^2);

    % Axil direction vectors
    % nde_1 = [delta_x; delta_y; delta_z]/Le; % axil 1
    vd1 = (n2 - n1)'/norm(n2 - n1); % eje Xele vector columna
    vd3 = cross(vd1, vaux)/norm(cross(vd1, vaux)); % lateral 3
    vd2 = cross(vd3, vd1); % lateral 2
    vdire = [vd1, vd2, vd3]; % N_e{e} = Ndir_e;  % almacen de vectores directores
    
    % DENTRO DE Geom_viga3D
    % Re matriz de cambio de base global-local corregida
    CXx = vd1'*dirX; CYx = vd1'*dirY; CZx = vd1'*dirZ; % Proyección x-local
    CXy = vd2'*dirX; CYy = vd2'*dirY; CZy = vd2'*dirZ; % Proyección y-local
    CXz = vd3'*dirX; CYz = vd3'*dirY; CZz = vd3'*dirZ; % Proyección z-local
    
    Re = [CXx CYx CZx; 
          CXy CYy CZy; 
          CXz CYz CZz];

    % Te matriz de transformaciones global-local del elemento
    R0 = zeros (3);
    Te = [Re R0 R0 R0; R0 Re R0 R0; R0 R0 Re R0; R0 R0 R0 Re];

end

% 3D beam stiffness (rigidez)
function [Keviga] = Keviga3D(E,G,Le,A,J,Iy,Iz)
    % rigideces escalares
    k1ax = E*A/Le ; k1tor = G*J/Le;
    k2fz = 12*E*Iz/Le^3; kf3z = 6*E*Iz/Le^2; kf4z = 4*E*Iz/Le; kf5z = 2*E*Iz/Le;
    kf2y = 12*E*Iy/Le^3; kf3y = 6*E*Iy/Le^2; kf4y = 4*E*Iy/Le; kf5y = 2*E*Iy/Le;
    % rigidez matricial a través de las submatrices nodales
    K11 = [k1ax 0 0 0 0 0; 0 k2fz 0 0 0 kf3z; 0 0 kf2y 0 -kf3y 0; 0 0 0 k1tor 0 0; 0 0 -kf3y 0 kf4y 0; 0 kf3z 0 0 0 kf4z];
    K12 = [-k1ax 0 0 0 0 0; 0 -k2fz 0 0 0 kf3z; 0 0 -kf2y 0 -kf3y 0; 0 0 0 -k1tor 0 0; 0 0 kf3y 0 kf5y 0; 0 -kf3z 0 0 0 kf5z];
    K22 = [k1ax 0 0 0 0 0; 0 k2fz 0 0 0 -kf3z; 0 0 kf2y 0 kf3y 0; 0 0 0 k1tor 0 0; 0 0 kf3y 0 kf4y 0; 0 -kf3z 0 0 0 kf4z];

    % MATRIZ DE RIGIDEZ LOCAL 3D
    Keviga = [K11 K12; K12' K22];

    
end

%% REPRESENTACIÓN GRÁFICA DEL SKYLINE (OCUPACIÓN DE K)
% Definimos el número total de grados de libertad basándonos en tu matriz K
GDL = size(K, 1); 
verK = zeros(GDL);

for s = 1:size(K, 1)
    for t = 1:size(K, 1)
        if abs(K(s, t)) > 1e-10
            verK(s, t) = 1; % Marcamos como 1 si hay rigidez (punto negro)
        else
            verK(s, t) = 0; % Marcamos como 0 si es vacío (punto blanco)
        end
    end
end

% Generación del gráfico
figure('Name', 'Skyline de la Matriz de Rigidez');
imagesc(verK);            % Crea el mapa de bits
colormap(flipud(gray));  % Invalida los colores para que sea Blanco (0) y Negro (1)
title('Ocupación de la Matriz de Rigidez K', 'FontSize', 12);
xlabel('Grados de Libertad (Columnas)', 'FontSize', 10);
ylabel('Grados de Libertad (Filas)', 'FontSize', 10);
axis square;             % Fuerza a que la matriz se vea cuadrada
grid on;                 % Ayuda a ver las divisiones de los nodos



%% REPRESENTACIÓN GRÁFICA DE LA ESTRUCTURA Y DEFORMADA
% Factor de escala para amplificar los desplazamientos
escala = 30; 

% Creamos una matriz con las coordenadas de todos los nodos para facilitar el acceso
n_coords = [n1; n2; n3; n4; n5; n6; n7];
% Matriz de conectividad (qué nodos forman cada uno de los 6 elementos)
elementos_conect = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7];

figure('Color', 'w', 'Name', 'Deformada Elástica vs Estructura Original');
hold on; grid on; axis equal; view(135, 30);

for e = 1:6
    % 1. Obtener índices de los nodos del elemento
    idx1 = elementos_conect(e,1); 
    idx2 = elementos_conect(e,2);
    
    % 2. Coordenadas originales de los nodos
    orig_n1 = n_coords(idx1, :);
    orig_n2 = n_coords(idx2, :);
    
    % Dibujar estructura original (línea negra discontinua)
    plot3([orig_n1(1) orig_n2(1)], [orig_n1(2) orig_n2(2)], [orig_n1(3) orig_n2(3)], 'k--');
    
    % 3. Extraer desplazamientos (Ux, Uy, Uz) de los dos nodos del elemento
    % GDLs Nodo 1: (idx-1)*6 + [1,2,3]
    u_n1 = U((idx1-1)*6 + (1:3))';
    u_n2 = U((idx2-1)*6 + (1:3))';
    
    % 4. Calcular coordenadas deformadas
    def_n1 = orig_n1 + escala * u_n1;
    def_n2 = orig_n2 + escala * u_n2;
    
    % Dibujar estructura deformada (línea azul con círculos)
    plot3([def_n1(1) def_n2(1)], [def_n1(2) def_n2(2)], [def_n1(3) def_n2(3)], 'b-o', 'LineWidth', 2);
end

title(['Deformada Elástica (Escala x', num2str(escala), ')']);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Estructura Original', 'Estructura Deformada');


% Matriz de rigidez global (consultar página 37 de las diapositias)

K_quest = K

% vector de fuerzas nodales equivalentes (consultar página 30 de las diapositias)

Fneq_quest = Fneq'

% desplazamientos

U_quest = U

% fuerzas de reacción

% Reacciones

FR1_quest = R1
FR7_quest = R7

% Fuerza en el desplazamiento impuesto
F_uload = K(20, :) * U - Fneq(20);

% Fuerzas internas locales (Elemento 1, 4 y 6)
% Elemento 1 (TIENE carga distribuida)
Fe1_quest = (Ke_local1 * (Te1 * U(gdle1)) - fneqe1)';

% Elemento 4 (NO tiene carga distribuida)
Fe4_quest = (Ke_local4 * Te4 * U(gdle4))';

% Elemento 6 (NO tiene carga distribuida)
Fe6_quest = (Ke_local6 * Te6 * U(gdle6))';


% Fin del código