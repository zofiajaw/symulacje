import numpy as np

phi_0 = 0.1 
dp = 0.1
mc_steps = 100


N = 100
x0 = -3
xN = 3
h = (xN-x0)/N


x = [None]*N
for i in range(N):
    x[i]=x0+(i)*h

#Zadanie poczatkowej wartosci funkcji falowej 
phi = [None]*N
for i in range(N):
    phi[i]=phi_0
    
#Energia poczatkowa kazdej funkcji falowej oddzielnie
v=[None]*(N)
t=[None]*(N)

for i in range(1, N+1):
    v[i] = (0.5 * np.power(x[i],2)) 
    t[i] = (0.5 * phi[i]*(2*phi[i] - phi[i-1]-phi[i+1]))/np.power(h,2)
    
#szczegolne przypadki
v[0]= 0.5 * np.power(x[0],2)
t[0] = (0.5 * phi[0] * (2*phi[0]))/np.power(h,2)
v[N-1]=0.5 * np.power(x[N-1],2) 
t[N-1]= (0.5 * phi[N-1]*(2*phi[N-1] - phi[N-1]))/np.power(h,2)

#Energia poczatkowa calkowita
U=0
T=0
Phi=0
for i in range(N):
    U = U+v[i]*np.power(phi[i],2)
    T = T+t[i]
    Phi = Phi+ np.power(phi[i],2)
E=(U + T)/Phi

#petla po krokach Monte Carlo
for j in range(mc_steps):
    for i in range(N):
        #Losowanie indeksu
        indeks = int(np.random.rand()*(N))
        phi_old = phi[indeks]
        #losowanie nowej funkcji falowej
        phi_new = phi_old +(np.random.rand() - 0.5) * dp
        #obliczanie nowej energii
        d_phi = phi_new - phi_old
        d_phi_pw = np.power(phi_new,2)-np.power(phi_old,2)
        if (indeks == 0):
            dT = (d_phi_pw - d_phi * phi[indeks+1])/np.power(h,2)
            dU = d_phi_pw * v[indeks]  
        elif (indeks == N-1):
            dT = (d_phi_pw - d_phi * (phi[indeks-1]))/np.power(h,2)
            dU = d_phi_pw * v[indeks] 
        else:
            dT = (d_phi_pw - d_phi * (phi[indeks+1] + phi[indeks-1]))/np.power(h,2)
            dU = d_phi_pw * v[indeks]
        E_new = (T+U+dT+dU)/(Phi+d_phi_pw)
        print (E_new)
        #sprawdzenie warunku
        if (E_new < E):
            phi[indeks] = phi_new
            E = E_new
            U = dU + U
            T = dT + T
            
        i=i+1
    j=j+i

#Zapis do pliku
arr=[]
for i in range(N):
    b=phi[i]
    a=x[i]
    arr.append(a)
    arr.append(b)   
file = '' #ścieżka do pliku
with open(file, 'w') as plik:
    for i in range(0,2*N,2):
       # plik.write('x: ')
        plik.write(str(arr[i]))
       # plik.write('\t\t')
        plik.write('phi(x): ')
        plik.write(str(arr[i+1]))
        plik.write('\n')
plik.close()