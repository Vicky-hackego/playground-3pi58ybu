
import matplotlib.pyplot as plt
import numpy as np

max_t = 360 #durée de la simulation
sub = 1 #subdvision
integrator = 'RK4'#Methode d'intégration : Euler ou RK4
R0 = 4.0 #R0 : taux de reproduction initial
Rc = 0.8 #Rc : taux de reproduction  en confinement
Rd = 1.2 #Rd : taux de reproduction  en déconfinement
tc = 20 #jour du confinement
td = 60 #jour du déconfinement
te = 6 #temps d'incubation
ti = 6 #temps d'infection
th = 20 #temps d'hospitalisation
tr = 180 #durée d'immunisation
k = 0.1 #proportion détectée
dr0 = 0.05*k #taux de mortalité
ir = 0.001 #portion contaminée au départ 

# Part réelle d'hospitalisation des malades (Belgique : 50% des cas détectés)
hr = 0.5*k

# Seuil de dépassement des hopitaux (Belgique : 50% de 5.6 lits/1000hab occupés)
threshold= 0.5*0.0056/hr

t = np.linspace(1., max_t,max_t*sub)

def seihrd(v,time):
    s,e,i,h,r,d = v
    Rt = R0 if time<tc else Rc if time<td else Rd
    dr = dr0 if h<threshold else dr0*(1+4*(h-threshold))
    return [-Rt/ti*s*i+r/tr,-e/te+Rt/ti*s*i, e/te-i/ti,i/ti-h/th,h/th*(1-dr)-r/tr,h/th*dr]
    
def euler(f,v0,t):
    v = np.array(v0)
    res=np.expand_dims(np.array(v0),axis=0)
    dt=t[1]-t[0]
    for s in t[1:]:
        v=v+np.array(f(v,s))*dt
        res = np.concatenate((res, np.expand_dims(v, axis=0)), axis=0)
    return res
    
def rk4(f,v0,t):
    v = np.array(v0)
    res=np.expand_dims(np.array(v0),axis=0)
    dt=t[1]-t[0]
    for s in t[1:]:
        k1=np.array(f(v,s))*dt
        k2=np.array(f(v+k1/2,s+dt/2))*dt
        k3=np.array(f(v+k2/2,s+dt/2))*dt
        k4=np.array(f(v+k3,s+dt))*dt
        v=v+(k1+2*k2+2*k3+k4)/6
        res = np.concatenate((res, np.expand_dims(v, axis=0)), axis=0)
    return res

if integrator == 'Euler':
    s,e,i,h,r,d=euler(seihrd, (1-ir,ir,0,0,0,0), t).T
else:
    s,e,i,h,r,d=rk4(seihrd, (1-ir,ir,0,0,0,0), t).T

fig=plt.figure()
plt.subplot(211)
plt.plot(t,s,label='Susceptibles')
plt.plot(t,r,label='Rétablis')
plt.plot(t,d,label='Décès')
plt.legend(loc=0)
plt.grid()
plt.subplot(212)
plt.plot(t,e,label='Incubations')
plt.plot(t,i,label='Infectieux')
plt.plot(t,h,label='Hospitalisations')
plt.legend(loc=0)
plt.grid()
plt.show()
fig.savefig('output.png', dpi=fig.dpi)
for k,_ in enumerate(t):
    print("%d %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f"%(k,s[k],e[k],i[k],h[k],r[k],d[k]))
print("TECHIO> open -s /project/target index.html")
