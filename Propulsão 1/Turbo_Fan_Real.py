import matplotlib.pyplot as plt
import numpy as np

class Turbofan_real:
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
    
    def Tau_de_f(self):
        tau_f = pi_f**((gamma_c-1)/(gamma_c*e_f))
        return tau_f
    
    def N_de_f(self):
        n_f = (((pi_f)**((gamma_c-1)/(gamma_c)))-1)/(self.Tau_de_f()-1)
        return n_f
    
    def F(self):
        f = (self.Tau_de_lambida() - (self.Tau_de_r()*self.Tau_de_c()))/(((n_b*h_pr)/(c_pc*t_0))-self.Tau_de_lambida())
        return f
    
    def Tau_de_t(self):
        tau_t = 1 - (((1)/(n_m*(1+self.F())))*((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c() - 1 + (alpha*(self.Tau_de_f()-1))))
        return tau_t
    
    def Pi_de_t(self):
        pi_t = self.Tau_de_t()**((gamma_t)/((gamma_t-1)*e_t))
        return pi_t
    
    def N_de_t(self):
        n_t = (1-self.Tau_de_t())/(1-(self.Tau_de_t()**(1/e_t)))
        return n_t
    
    def P_t9IP_9(self):
        p_t9Ip_9 = p_0Ip_9*self.Pi_de_r()*self.Pi_de_d()*pi_c*pi_b*self.Pi_de_t()*pi_n
        return p_t9Ip_9 
    
    def M_9(self):
        m_9 = (((2)/(gamma_t-1))*(((self.P_t9IP_9())**((gamma_t-1)/(gamma_t)))-1))**(0.5)
        return m_9
    
    def T_9IT_0(self):
        t_9It_0 = ((self.Tau_de_t()*self.Tau_de_lambida())/(self.P_t9IP_9()**((gamma_t-1)/(gamma_t))))*((c_pc)/(c_pt))
        return t_9It_0
    
    def V_9Ia_0(self):
        v_9Ia_0 = self.M_9()*((((gamma_t*self.R_t())/(gamma_c*self.R_c()))*(self.T_9IT_0()))**(0.5))
        return v_9Ia_0
    
    def P_t19IP_19(self):
        p_t19Ip_19 = p_0Ip_19*self.Pi_de_r()*self.Pi_de_d()*pi_f*pi_fn
        return p_t19Ip_19

    def M_19(self):
        m_19 = (((2)/(gamma_c-1))*(((self.P_t19IP_19())**((gamma_c-1)/(gamma_c)))-1))**(0.5)
        return m_19 
    
    def T_19IT_0(self):
        t_19It_0 = (self.Tau_de_f()*self.Tau_de_r())/(self.P_t19IP_19()**((gamma_c-1)/(gamma_c)))
        return t_19It_0
    
    def V_19Ia_0(self):
        v_19Ia_0 = self.M_19()*((self.T_19IT_0())**(0.5))
        return v_19Ia_0
    
    def Empuxo(self):
        parte_1 = ((alpha)/(1+alpha))*(self.A_0())*((self.V_19Ia_0())-(m_0)+(((self.T_19IT_0())/(self.V_19Ia_0()))*((1-p_0Ip_19)/(gamma_c))))
        parte_2 = ((((1)/(1+alpha))*(self.A_0()))*(((1+self.F())*(self.V_9Ia_0()))-(m_0)+((1+self.F())*((self.R_t())/(self.R_c()))*((self.T_9IT_0())/(self.V_9Ia_0()))*((1-p_0Ip_9)/(gamma_c)))))
        empuxo = parte_1+parte_2
        return empuxo
    
    def S(self):
       s = self.F()/((1+alpha)*self.Empuxo())
       return s
    
    def FR(self):
       parte_1 = ((self.V_19Ia_0())-(m_0)+(((self.T_19IT_0())/(self.V_19Ia_0()))*((1-p_0Ip_19)/(gamma_c))))
       parte_2 = (((1+self.F())*(self.V_9Ia_0()))-(m_0)+((1+self.F())*((self.R_t())/(self.R_c()))*((self.T_9IT_0())/(self.V_9Ia_0()))*((1-p_0Ip_9)/(gamma_c)))) 
       fr = (parte_2)/(parte_1)
       return fr
    
    def N_t(self):
       n_ter = (((self.A_0())**(2))*(((1+self.F())*((self.V_9Ia_0())**2))+(alpha*((self.V_19Ia_0())**2))-((1+alpha)*(m_0**2))))/(2*self.F()*h_pr)
       return n_ter
    
    def N_p(self):
       n_p = (2*m_0*((((1+self.F())*((self.V_9Ia_0())))+(alpha*((self.V_19Ia_0())))-((1+alpha)*(m_0)))))/((((1+self.F())*((self.V_9Ia_0())**2))+(alpha*((self.V_19Ia_0())**2))-((1+alpha)*(m_0**2))))
       return n_p
    
    def N_o(self):
       n_o = self.N_t()*self.N_p()
       return n_o


def Nivel_motor(nivel):
    global t_t4, pi_dmax, pi_b, pi_n, e_c, e_t, n_b, e_f
    if nivel == 0:        
        pi_dmax = 0.99
        e_c = 0.90
        e_f = 0.89
        pi_b = 0.96
        n_b = 0.99
        e_t = 0.89
        pi_n = 0.99  
        t_t4 = 3000
        if convert:
            t_t4 = t_t4/1.8 


    elif nivel == 1:        
        pi_dmax = 0.90  #ajustavel
        e_t = 0.80  #ajustavel
        pi_n = 0.95  #ajustavel

        t_t4 = 1110 
        pi_b = 0.90
        e_c = 0.80
        e_f = 0.78
        n_b = 0.85
        if convert:
            t_t4 = 2000 
   

    elif nivel == 2:
        pi_dmax = 0.95  #ajustavel
        e_t = 0.83  #ajustavel
        pi_n = 0.97  #ajustavel
     
        e_c = 0.84
        e_f = 0.82
        pi_b = 0.92
        n_b = 0.91
        t_t4 = 1390  
        if convert:
            t_t4 = 2500 
  

    elif nivel == 3:
        pi_dmax = 0.98  #ajustavel
        e_t = 0.87  #ajustavel
        pi_n = 0.98  #ajustavel

        t_t4 = 1780         
        pi_b = 0.94
        e_c = 0.88
        e_f = 0.86
        n_b = 0.98
        if convert:
            t_t4 = 3200 
     
    
    elif nivel == 4:
        pi_dmax = 0.995  #ajustavel
        e_t = 0.89  #ajustavel
        pi_n = 0.995  #ajustavel

        t_t4 = 2000         
        pi_b = 0.95
        e_c = 0.90
        e_f = 0.89
        n_b = 0.99
        if convert:
            t_t4 = 3600

#####################################################################################

motor = Turbofan_real()
nivel_motor = [0]

# Variaveis de entrada:
lista_alpha = [8]
lista_pi_c = [2,40]

# Constantes de entrada em SI:
t_0 = 216.7         #[K]
t_t4 = 1670       #[K]   
c_pc = 1004     #[J/kg K]
c_pt = 1096    #[J/kg K]
h_pr = 42800000     #[J/kg]

gamma_c = 1.4 
gamma_t = 1.33 
pi_f = 1.7
pi_fn = 0.99
n_m = 0.99
p_0Ip_9 = 0.9
p_0Ip_19 = 0.9
m_0 = 0.8
pi_c_ = 36
alpha = lista_alpha[0]

convert = True
if convert:
    # Constantes de entrada FORA do SI:
    t_0 = 390         #[R]
    t_t4 = 3000       #[R]   
    c_pc = 0.240     #[Btu/lbm R]
    c_pt = 0.276    #[Btu/lbm R]
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
    dado = format(pi_f, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_f",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_fn, ".3f")
    print("%s\t%s\t\t%s" % ("| pi_{fn}",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_b, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_b",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_n, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_n",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_f, ".3f")
    print("%s\t\t%s\t\t%s" % ("| e_f",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_c, ".3f")
    print("%s\t\t%s\t\t%s" % ("| e_c",dado," "))
    print("|-------------------------------------------|")
    dado = format(e_t, ".3f")
    print("%s\t\t%s\t\t%s" % ("| e_t",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_b, ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_b",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_m, ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_m",dado," "))
    print("|-------------------------------------------|")
    dado = format(p_0Ip_9, ".3f")
    print("%s\t%s\t\t%s" % ("| p_0/p_9",dado," "))
    print("|-------------------------------------------|")
    dado = format(p_0Ip_19, ".3f")
    print("%s\t%s\t\t%s" % ("| p_0/p_19",dado," "))
    print("|-------------------------------------------|")
    dado = format(m_0, ".3f")
    print("%s\t\t%s\t\t%s" % ("| m_0",dado," "))
    print("|-------------------------------------------|")
    print(" ")


# Graficos:
plt.figure(figsize = ((15, 5)))
pi_c = np.linspace(min(lista_pi_c) ,max(lista_pi_c))  #Limites do eixo X

# Grafico 1 ----------------------------------------------------------
plt.subplot(2, 2, 1)

for i in range(len(nivel_motor)):
    Nivel_motor(nivel_motor[i])
    plt.plot(pi_c, motor.Empuxo(), label=r'$tec. level = ${}'.format(nivel_motor[i]))

            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(0, 800)
plt.xlabel(r'$\pi_c$' , fontsize=15)
plt.ylabel(r'$F/ \dot m$   $[N / (kg /s)]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 2 ----------------------------------------------------------
plt.subplot(2, 2, 2)
for i in range(len(nivel_motor)):
    Nivel_motor(nivel_motor[i])
    plt.plot(pi_c, motor.S()*1000000, label=r'$tec. level = ${}'.format(nivel_motor[i]))
            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(15, 50)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$S$   $[(mg/s) / N]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 3 ----------------------------------------------------------
plt.subplot(2, 2, 3)
for i in range(len(nivel_motor)):
    Nivel_motor(nivel_motor[i])
    plt.plot(pi_c, motor.F(), label=r'$tec. level = ${}'.format(nivel_motor[i]))
            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$f$' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 4 ----------------------------------------------------------


# Plot Grafico ----------------------------------------------------------
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
plt.show()


    
# Tabela de saida:

pi_c = pi_c_
for i in range(len(nivel_motor)):
    print(" ")
    Nivel_motor(nivel_motor[i])
    print(" \t\tTEC LEVEL = {}".format(nivel_motor[i]))
    print("|-------------------------------------------|")
    print("VALORES INICIAIS")
    dado = format(alpha, ".3f")
    print("%s\t\t%s\t\t%s" % ("| alpha",dado," "))
    print("|-------------------------------------------|")
    print("%s\t\t%s\t\t%s" % ("| pi_c",pi_c," "))
    print("|-------------------------------------------|")
    print(" ")
    print(" ")
    print(" SAIDA\t\tVALOR\t\tUNIDADE")
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
    dado = format(motor.FR(), ".3f")
    print("%s\t\t%s\t%s" % ("| FR",dado," "))
    print("|-------------------------------------------|")
    dado = format(motor.N_de_t(), ".3f")
    print("%s\t\t%s\t%s" % ("| eta_t",dado," "))
    print("|-------------------------------------------|")
    dado = format(motor.N_de_c(), ".3f")
    print("%s\t\t%s\t%s" % ("| eta_c",dado," "))
    print("|-------------------------------------------|")
    dado = format(motor.N_p(), ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_P",dado," "))
    print("|-------------------------------------------|")
    dado = format(motor.N_t(), ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_T",dado," "))
    print("|-------------------------------------------|")
    dado = format(motor.N_o(), ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_O",dado," "))
    print("|-------------------------------------------|")
    print(" ")
    print(" ")




