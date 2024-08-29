import matplotlib.pyplot as plt
import numpy as np

class Turbofan_prop_real:
    def __init__(self):
        s = 1

    def R_c(self):
        r_c = ((gamma_c-1)/(gamma_c))*(c_pc)
        return r_c
    
    def R_t(self):
        r_t = ((gamma_t-1)/(gamma_t))*(c_pt)
        return r_t

    def A_0(self):
        a_0 = (gamma_c*self.R_c()*t_0)**(0.5)
        return a_0
    
    def V_0(self):
        v_0 = self.A_0()*m_0
        return v_0
    
    def Tau_de_r(self):
        tau_r = 1 + (((gamma_c - 1)/(2))*(m_0**2))
        return tau_r
    
    def Pi_de_r(self):
        pi_r = (self.Tau_de_r())**((gamma_c)/(gamma_c-1))
        return pi_r
    
    def Pi_de_d(self):
        if m_0 <= 1:
            n_r = 1
        else:
            n_r = 1 - (0.075*((m_0-1)**1.35))
        pi_d = pi_dmax*n_r
        return pi_d
    
    def Tau_de_lambida(self):
        tau_lambida = (t_t4*c_pt)/(t_0*c_pc)
        return tau_lambida
    
    def Tau_de_c(self):
        tau_c = pi_c**((gamma_c-1)/(gamma_c*e_c))
        return tau_c
    
    def N_de_c(self):
        n_c = (((pi_c)**((gamma_c-1)/(gamma_c)))-1)/(self.Tau_de_c()-1)
        return n_c
    
    def F(self):
        f = (self.Tau_de_lambida() - (self.Tau_de_r()*self.Tau_de_c()))/(((n_b*h_pr)/(c_pc*t_0))-self.Tau_de_lambida())
        return f
    
    def Tau_de_th(self):
        tau_th = 1-(((1)/(n_mh*(1+self.F())))*((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1))
        return tau_th
    
    def Pi_de_th(self):
        pi_th = (self.Tau_de_th())**((gamma_t)/((gamma_t-1)*(e_th)))
        return pi_th
    
    def N_de_th(self):
        n_th = (1-self.Tau_de_th())/(1-((self.Tau_de_th())**((1)/(e_th))))
        return n_th
    
    def A(self):
        a = (((gamma_c-1)/(2))*((m_0**2)/(self.Tau_de_lambida()*self.Tau_de_th())))/((n_prop*n_g*n_ml)**(2))
        return a
    
    def II(self):
        ii = (self.Pi_de_r()*self.Pi_de_d()*pi_c*pi_b*pi_n)**((gamma_t-1)/(gamma_t))
        return ii
    
    def Tau_de_tl(self):
        if com_tau_t:
            tau_tl = (tau_t)/(self.Tau_de_th())
        else:
            trava = True
            tau_tl = (((self.Tau_de_th())**((-1)/(e_th)))/(self.II())) + self.A()
            tau_t_antigo = tau_tl
            while trava:
                tau_tl = ((((self.Tau_de_th())**((-1)/(e_th)))/(self.II()))*((tau_t_antigo)**((e_tl-1)/(e_tl))))+((self.A())*((1+(((1-e_tl)/(e_tl))*((((self.Tau_de_th())**((-1)/(e_th)))*((tau_t_antigo)**((-1)/(e_tl))))/(self.II()))))**(2)))
                teste = ((tau_tl - tau_t_antigo)**(2))**(0.5)
                if np.all(teste <= 0.0001):
                    trava = False
                tau_t_antigo = tau_tl
        return tau_tl
    
    
    def Pi_de_tl(self):
        pi_tl = (self.Tau_de_tl())**((gamma_t)/((gamma_t-1)*e_tl))
        return pi_tl
    
    def N_de_tl(self):
        n_tl = (1-self.Tau_de_tl())/(1-((self.Tau_de_tl())**((1)/(e_tl))))
        return n_tl
    
    def P_t9IP_0(self):
        p_t9Ip_0 =(self.Pi_de_r()*self.Pi_de_d()*pi_c*pi_b*self.Pi_de_th()*self.Pi_de_tl()*pi_n)
        return p_t9Ip_0 
    
    def M_9(self):
        parte_1 = self.P_t9IP_0()
        parte_2 = ((gamma_t+1)/(2))**((gamma_t)/(gamma_t-1))
        if np.all(parte_1 > parte_2):
            m_9 = 1
        else:
            m_9 = (((2)/(gamma_t-1))*(((self.P_t9IP_0)**((gamma_t-1)/(gamma_t)))-1))**(0.5)
        return m_9
    
    def P_t9IP_9(self):
        parte_1 = self.P_t9IP_0()
        parte_2 = ((gamma_t+1)/(2))**((gamma_t)/(gamma_t-1))
        if np.all(parte_1 > parte_2):
            p_t9Ip_9 = parte_2
        else:
            p_t9Ip_9 = self.P_t9IP_0()
        return p_t9Ip_9
    
    def P_0IP_9(self):
        parte_1 = self.P_t9IP_0()
        parte_2 = ((gamma_t+1)/(2))**((gamma_t)/(gamma_t-1))
        if np.all(parte_1 > parte_2):
            p_0Ip_9 =(self.P_t9IP_9())/(self.P_t9IP_0())
        else:
            p_0Ip_9 = 1
        return p_0Ip_9 
    
    def T_t9IT_0(self):
        t_t9It_0 = ((t_t4)/(t_0))*(self.Tau_de_th())*(self.Tau_de_tl())
        return t_t9It_0
    
    def T_9IT_0(self):
        t_9It_0 = (self.T_t9IT_0())/((self.P_t9IP_9())**((gamma_t-1)/(gamma_t)))
        return t_9It_0
    
    def V_9Ia_0(self):
        v_9Ia_0 = (((2*self.Tau_de_lambida()*self.Tau_de_th()*self.Tau_de_tl())/(gamma_c-1))*(1-((self.P_t9IP_9())**((1-gamma_t)/(gamma_t)))))**(0.5)
        return v_9Ia_0
    
    def C_prop(self):
        c_prop = n_prop*n_g*n_ml*(1+self.F())*self.Tau_de_lambida()*self.Tau_de_th()*(1-self.Tau_de_tl())
        return c_prop
    
    def C_c(self):
        c_c = (gamma_c-1)*m_0*(((1+self.F())*(self.V_9Ia_0()))-m_0+((1+self.F())*((self.R_t())/(self.R_c()))*((self.T_9IT_0())/(self.V_9Ia_0()))*((1-self.P_0IP_9())/(gamma_c))))
        return c_c
    
    def C_tot(self):
        c_tot = self.C_prop() + self.C_c()
        return c_tot
    
    def Empuxo(self):
        empuxo = (self.C_tot()*c_pc*t_0)/(self.V_0())
        return empuxo
    
    def S(self):
       s = self.F()/(self.Empuxo())
       return s
    
    def W_m(self):
       w_m = self.C_tot()*c_pc*t_0
       return w_m 
    
    def S_p(self):
       s_p = (self.F())/(self.C_tot()*c_pc*t_0)
       return s_p
    
    def N_t(self):
       n_ter = (self.C_tot())/((self.F()*h_pr)/(c_pc*t_0))
       return n_ter
    
    def N_p(self):
       n_p = (self.C_tot())/(((self.C_prop())/(n_prop))+(((gamma_c-1)/(2))*(((1+self.F())*((self.V_9Ia_0())**2))-(m_0**2))))
       return n_p
    
def Nivel_motor(nivel):
    global t_t4, pi_dmax, pi_b, pi_n, e_c, e_t, n_b
    if nivel == 0:        
        pi_dmax = 0.98
        e_c = 0.92 
        pi_b = 0.98
        n_b = 0.99
        e_t = 0.91 
        pi_n = 0.98  
        t_t4 = 1800
        if convert:
            t_t4 = t_t4/1.8 


    elif nivel == 1:        
        pi_dmax = 0.90  #ajustavel
        e_t = 0.80  #ajustavel
        pi_n = 0.95  #ajustavel

        t_t4 = 1110 
        pi_b = 0.90
        e_c = 0.80
        n_b = 0.85
        if convert:
            t_t4 = 2000 
   

    elif nivel == 2:
        pi_dmax = 0.95  #ajustavel
        e_t = 0.85  #ajustavel
        pi_n = 0.97  #ajustavel

        t_t4 = 1390       
        pi_b = 0.92
        e_c = 0.84
        n_b = 0.91
        if convert:
            t_t4 = 2500 
  

    elif nivel == 3:
        pi_dmax = 0.94  #ajustavel
        e_t = 0.87  #ajustavel
        pi_n = 0.95  #ajustavel

        t_t4 = 1780         
        pi_b = 0.94
        e_c = 0.88
        n_b = 0.98
        if convert:
            t_t4 = 3200 
     
    
    elif nivel == 4:
        pi_dmax = 0.995  #ajustavel
        e_t = 0.90  #ajustavel
        pi_n = 0.995  #ajustavel

        t_t4 = 2000         
        pi_b = 0.95
        e_c = 0.90
        n_b = 0.99
        if convert:
            t_t4 = 3600   
    
#####################################################################################

motor = Turbofan_prop_real()
nivel_motor = [2]

# Variaveis de entrada:
lista_pi_c = [0.001,14]
lista_tau_t = [0.5]

# Constantes de entrada em SI:
t_0 = 250         #[K]
c_pc = 1004     #[J/kg K]
c_pt = 1108    #[J/kg K]
h_pr = 42800000     #[J/kg]

gamma_c = 1.4 
gamma_t = 1.35 
e_tl = 0.91
e_th = 0.89
n_prop = 0.83 
n_g = 0.99 
n_mh = 0.99
n_ml = 0.99
m_0 = 0.7 
pi_c_ = 4.6

com_tau_t = False

convert = True
if convert:
    # Constantes de entrada FORA do SI:
    t_0 = 447         #[R]
    t_t4 = 3000       #[R]   
    c_pc = 0.240     #[Btu/lbm R]
    c_pt = 0.265    #[Btu/lbm R]
    h_pr = 18400     #[Btu/lbm]
    
    #conversão
    t_0 = t_0/1.8             #[K]
    t_t4 = t_t4/1.8           #[K]
    c_pc = c_pc  * 4.1868   #[kJ/kg K]
    c_pc = c_pc * 1000        #[J/kg K]
    c_pt = c_pt  * 4.1868   #[kJ/kg K]
    c_pt = c_pt * 1000        #[J/kg K]
    h_pr = h_pr * 2.326       #[kJ/kg]
    h_pr = h_pr * 1000        #[J/kg]


# Tabela com os dados de entrada
for i in range(len(nivel_motor)):
    print(" ")
    Nivel_motor(nivel_motor[i])
    print(" \t\tTEC LEVEL = {}".format(nivel_motor[i]))
    print("|-------------------------------------------|")
    print(" ENTRADA\tVALOR\t\tUNIDADE")
    print("|-------------------------------------------|")
    dado = format(t_0, ".3f")
    print("%s\t\t%s\t\t%s" % ("| T_0",dado,"K"))
    print("|-------------------------------------------|")
    dado = format(t_t4, ".3f")
    print("%s\t%s\t%s" % ("| T_{t4}",dado,"K"))
    print("|-------------------------------------------|")
    dado = format(c_pc, ".3f")
    print("%s\t%s\t%s" % ("| c_{pc}",dado,"[J/kg K]"))
    print("|-------------------------------------------|")
    dado = format(c_pt, ".3f")
    print("%s\t%s\t%s" % ("| c_{pt}",dado,"[J/kg K]"))
    print("|-------------------------------------------|")
    dado = format(h_pr, ".3f")
    print("%s\t%s\t%s" % ("| h_{pr}",dado,"[J/kg]"))
    print("|-------------------------------------------|")
    dado = format(gamma_c, ".3f")
    print("%s\t%s\t\t%s" % ("| gamma_c",dado," "))
    print("|-------------------------------------------|")
    dado = format(gamma_t, ".3f")
    print("%s\t%s\t\t%s" % ("| gamma_t",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_dmax, ".3f")
    print("%s\t%s\t\t%s" % ("| pi_{dmax}",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_b, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_b",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_n, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_n",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_th, ".3f")
    print("%s\t%s\t\t%s" % ("| e_{tH}",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_c, ".3f")
    print("%s\t\t%s\t\t%s" % ("| e_c",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_tl, ".3f")
    print("%s\t%s\t\t%s" % ("| e_{tL}",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_b, ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_b",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_mh, ".3f")
    print("%s\t%s\t\t%s" % ("| eta_{mH}",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_ml, ".3f")
    print("%s\t%s\t\t%s" % ("| eta_{mL}",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_prop, ".3f")
    print("%s\t%s\t\t%s" % ("| eta_{prop}",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_g, ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_g",dado," "))
    print("|-------------------------------------------|")
    dado = format(m_0, ".3f")
    print("%s\t\t%s\t\t%s" % ("| M_0",dado," "))
    print("|-------------------------------------------|")
    print(" ")
    print("ººººººººººººººººººººººººººººººººººººººººººººº")


# Graficos:
plt.figure(figsize = ((15, 5)))
pi_c = np.linspace(min(lista_pi_c) ,max(lista_pi_c))  #Limites do eixo X

# Grafico 1 ----------------------------------------------------------
plt.subplot(1, 2, 1)

if com_tau_t:
    for i in range(len(lista_tau_t)):
        tau_t = lista_tau_t[i]
        plt.plot(pi_c, motor.Empuxo(), label=r'$\tau_t = ${}'.format(tau_t))
else:
    for i in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.Empuxo(), label=r'$tec. level = ${}'.format(nivel_motor[i]))
            
y_ = 120*9.8067
plt.axhline(y=y_, color='red', linestyle='--')
x_ = 1.5
plt.axvline(x=x_, color='red', linestyle='--')
plt.scatter(x_, y_, color='red')
plt.text(x_, y_, f'({x_:.1f}, {y_:.1f})', fontsize=12, ha='left', va='bottom')
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(0, 3700)
plt.xlabel(r'$\pi_c$' , fontsize=15)
plt.ylabel(r'$F/ \dot m$   $[N / (kg /s)]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 2 ----------------------------------------------------------
plt.subplot(1, 2, 2)
if com_tau_t:
    for i in range(len(lista_tau_t)):
        tau_t = lista_tau_t[i]
        plt.plot(pi_c, motor.S()*1000000, label=r'$\tau_t = ${}'.format(tau_t))
else:
    for i in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.S()*1000000, label=r'$tec. level = ${}'.format(nivel_motor[i]))    

y_ = 0.8*28.325
plt.axhline(y=y_, color='red', linestyle='--')
x_ = 4.6
plt.axvline(x=x_, color='red', linestyle='--')
plt.scatter(x_, y_, color='red')
plt.text(x_, y_, f'({x_:.1f}, {y_:.1f})', fontsize=12, ha='left', va='bottom')      
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(15, 40)
plt.xlim(2, 14)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$S$   $[(mg/s) / N]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Plot Grafico ----------------------------------------------------------
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
plt.show()


    
# Tabela de saida:
for j in range(len(lista_tau_t)):
    if com_tau_t:
        tau_t = lista_tau_t[j]
        printa_tau_T = False
    else:
        j = len(lista_tau_t)
        printa_tau_T = True
    pi_c = pi_c_
    for i in range(len(nivel_motor)):
        print(" ")
        Nivel_motor(nivel_motor[i])
        print(" \t\tTEC LEVEL = {}".format(nivel_motor[i]))
        print("|-------------------------------------------|")
        print("VALORES INICIAIS")
        if com_tau_t:
            dado = format(tau_t , ".3f")
            print("%s\t%s\t\t%s" % ("| tau_t ",dado," "))
            print("|-------------------------------------------|")
        print("%s\t\t%s\t\t%s" % ("| pi_c",pi_c," "))
        print("|-------------------------------------------|")
        print(" ")
        print(" ")
        print(" SAIDA\t\tVALOR\t\tUNIDADE")
        print("|-------------------------------------------|")
        if printa_tau_T:
            dado = format(motor.Tau_de_tl() , ".3f")
            print("%s\t\t%s\t\t%s" % ("| tau_tl* ",dado," "))
            print("|-------------------------------------------|")
        dado = format(motor.S()*1000000, ".3f")
        print("%s\t\t%s\t\t%s" % ("| S",dado,"mg s/N"))
        print("|-------------------------------------------|")
        dado = format(motor.Empuxo(), ".3f")
        if motor.Empuxo() <1000:
            print("%s\t\t%s\t\t%s" % ("| F/m_0",dado,"N s/kg"))
            print("|-------------------------------------------|")
        else:
            print("%s\t\t%s\t%s" % ("| F/m_0",dado,"N s/kg"))
            print("|-------------------------------------------|")
        dado = format(motor.F(), ".3f")
        print("%s\t\t%s\t%s" % ("| f",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.W_m(), ".3f")
        print("%s\t\t%s\t%s" % ("| W/m_0",dado,"W s/kg  "))
        print("|-------------------------------------------|")
        dado = format(motor.S_p(), ".3f")
        print("%s\t\t%s\t\t%s" % ("| S_p",dado,"mg /W s "))
        print("|-------------------------------------------|")
        dado = format(motor.C_c(), ".3f")
        print("%s\t\t%s\t%s" % ("| C_c",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.C_prop(), ".3f")
        print("%s\t%s\t%s" % ("| C_{prop}",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.C_tot(), ".3f")
        print("%s\t%s\t%s" % ("| C_{tot}",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.N_p(), ".3f")
        print("%s\t\t%s\t\t%s" % ("| eta_P",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.N_t(), ".3f")
        print("%s\t\t%s\t\t%s" % ("| eta_T",dado," "))
        print("|-------------------------------------------|")
        print(" ")
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")








