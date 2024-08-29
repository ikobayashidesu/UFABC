import matplotlib.pyplot as plt
import numpy as np

class Turbofan_mist_real:
    def __init__(self):
        s = 1

    def R_c(self):
        r_c = ((gamma_c-1)/(gamma_c))*(c_pc)
        return r_c
    
    def R_t(self):
        r_t = ((gamma_t-1)/(gamma_t))*(c_pt)
        return r_t

    def R_ab(self):
        r_ab = ((gamma_ab-1)/(gamma_ab))*(c_pab)
        return r_ab

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
    
    def Tau_de_lambida_AB(self):
        tau_lambida_ab = (t_t7*c_pab)/(t_0*c_pc)
        return tau_lambida_ab
    
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
    
    def Alpha(self):
        alpha = (((n_m)*(1+self.F())*((self.Tau_de_lambida())/(self.Tau_de_r()))*(1-(((pi_f)/(pi_c*pi_b))**((e_t*(gamma_t-1))/(gamma_t)))))-(self.Tau_de_c())+(1))/(self.Tau_de_f()-1)
        return alpha
    
    def Tau_de_t(self):
        tau_t = 1 - (((1)/(n_m*(1+self.F())))*((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c() - 1 + (self.Alpha()*(self.Tau_de_f()-1))))
        return tau_t
    
    def Pi_de_t(self):
        pi_t = self.Tau_de_t()**((gamma_t)/((gamma_t-1)*e_t))
        return pi_t
    
    def N_de_t(self):
        n_t = (1-self.Tau_de_t())/(1-(self.Tau_de_t()**(1/e_t)))
        return n_t
    
    def P_t16IP_16(self):
        p_t16Ip_16 = (pi_f)/(pi_c*pi_b*self.Pi_de_t())
        return p_t16Ip_16
    
    def M_16(self):
        m_16 = (((2)/(gamma_t-1))*((((self.P_t16IP_16())*((1+(((gamma_t-1)/(2))*((m_6)**(2))))**((gamma_t)/(gamma_t-1))))**((gamma_c-1)/(gamma_c)))-1))**(0.5)
        return m_16
    
    def Alpha_l(self):
        alpha_l = (self.Alpha())/(1+self.F())
        return alpha_l
    
    def C_p6a(self):
        c_p6a = (c_pt+(self.Alpha_l()*c_pc))/(1+self.Alpha_l())
        return c_p6a
    
    def R_6a(self):
        r_6a = (self.R_t()+(self.Alpha_l()*self.R_c()))/(1+self.Alpha_l())
        return r_6a
    
    def Gamma_6a(self):
        gamma_6a = (self.C_p6a())/(self.C_p6a()-self.R_6a())
        return gamma_6a
    
    def T_t16IT_6(self):
        t_t16It_t6 = (t_0*self.Tau_de_r()*self.Tau_de_f())/(t_t4*self.Tau_de_t())
        return t_t16It_t6
    
    def Tau_de_M(self):
        tau_m = ((c_pt)/(self.C_p6a()))*((1+(self.Alpha_l()*((c_pc)/(c_pt))*(self.T_t16IT_6())))/(1+self.Alpha_l()))
        return tau_m
    
    def Fluxo_6(self):
        fluxo_6 = ((m_6**2)*(1+(((gamma_t-1)/(2))*(m_6**2))))/((1+(gamma_t*(m_6**2)))**2)
        return fluxo_6
    
    def Fluxo_16(self):
        fluxo_16 = ((self.M_16()**2)*(1+(((gamma_c-1)/(2))*(self.M_16()**2))))/((1+(gamma_c*(self.M_16()**2)))**2)
        return fluxo_16
    
    def Fluxo_6a(self):
        fluxo_6a = (((1+self.Alpha_l())/(((1)/(self.Fluxo_6()**(0.5)))+(self.Alpha_l()*(((self.R_c()/self.R_t())*((gamma_t)/(gamma_c))*((self.T_t16IT_6())/(self.Fluxo_16())))**(0.5)))))**(2))*((self.R_6a())/(self.R_t()))*((gamma_t)/(self.Gamma_6a()))*(self.Tau_de_M())
        return fluxo_6a
    
    def M_6a(self):
        M_6a = ((2*self.Fluxo_6a())/(1-(2*self.Gamma_6a()*self.Fluxo_6a())+((1-(2*(self.Gamma_6a()+1)*self.Fluxo_6a()))**(0.5))))**(0.5)
        return M_6a
    
    def A_16IA_6(self):
        a_16Ia_6 = (self.Alpha_l()*((self.T_t16IT_6())**(0.5)))/((self.M_16()/m_6)*((((gamma_c)/(gamma_t))*(self.R_t()/self.R_c())*((1+(((gamma_c-1)/(2))*(self.M_16()**2)))/(1+(((gamma_t-1)/(2))*(m_6**2)))))**(0.5)))
        return a_16Ia_6
    
    def Pi_de_m_ideal(self):
        mfp_6 = (m_6*(((gamma_t)/(self.R_t()))**(0.5)))/((1+(((gamma_t-1)/(2))*(m_6**2)))**((gamma_t+1)/(2*(gamma_t-1))))
        mfp_6a = (self.M_6a()*(((self.Gamma_6a())/(self.R_6a()))**(0.5)))/((1+(((self.Gamma_6a()-1)/(2))*(self.M_6a()**2)))**((self.Gamma_6a()+1)/(2*(self.Gamma_6a()-1))))
        pi_m_ideal = (((1+self.Alpha_l())*(self.Tau_de_M()**(0.5)))/(1+self.A_16IA_6()))*((mfp_6)/(mfp_6a))
        return pi_m_ideal
    
    def Pi_de_m(self):
        pi_m = pi_mmax*self.Pi_de_m_ideal()
        return pi_m
    
    def P_t9IP_9(self):
        p_t9Ip_9 = p_0Ip_9*self.Pi_de_r()*self.Pi_de_d()*pi_c*pi_b*self.Pi_de_t()*pi_n*self.Pi_de_m()*pi_ab
        return p_t9Ip_9 
    
    def C_p9(self):
        if after_burn:
            c_p9 = c_pab
        else:
            c_p9 = self.C_p6a()
        return c_p9
    
    def R_9(self):
        if after_burn:
            r_9 = self.R_ab()
        else:
            r_9 = self.R_6a()
        return r_9
    
    def Gamma_9(self):
        if after_burn:
            gamma_9 = gamma_ab
        else:
            gamma_9 = self.Gamma_6a()
        return gamma_9
    
    def F_AB(self):
        if after_burn:
            f_ab = (1+((self.F())/(1+self.Alpha())))*(((self.Tau_de_lambida_AB())-(((self.C_p6a())/(c_pt))*(self.Tau_de_lambida())*(self.Tau_de_t())*(self.Tau_de_M())))/(((n_ab*h_pr)/(c_pc*t_0))-(self.Tau_de_lambida_AB())))
        else:
            f_ab = 0
        return f_ab
    
    def T_9IT_0(self):
        if after_burn:
            t_9It_0 = ((t_t7)/(t_0))/((self.P_t9IP_9())**((self.Gamma_9()-1)/(self.Gamma_9()))) 
        else:
            t_9It_0 = ((t_t4*self.Tau_de_t()*self.Tau_de_M())/(t_0))/((self.P_t9IP_9())**((self.Gamma_9()-1)/(self.Gamma_9())))
        
        return t_9It_0
    
    def M_9(self):
        m_9 = (((2)/(self.Gamma_9()-1))*(((self.P_t9IP_9())**((self.Gamma_9()-1)/(self.Gamma_9())))-1))**(0.5)
        return m_9
    
    def V_9Ia_0(self):
        v_9Ia_0 = self.M_9()*((((self.Gamma_9()*self.R_9())/(gamma_c*self.R_c()))*(self.T_9IT_0()))**(0.5))
        return v_9Ia_0
    
    def F_o(self):
        f_o = ((self.F())/(1+self.Alpha()))+(self.F_AB())
        return f_o
    
    def Empuxo(self):
        empuxo = self.A_0()*(((1+self.F_o())*(self.V_9Ia_0()))-(m_0)+((1+self.F_o())*((self.R_9())/(self.R_c()))*(self.T_9IT_0()/self.V_9Ia_0())*((1-p_0Ip_9)/(gamma_c))))
        return empuxo
    
    def S(self):
       s = (self.F_o()/(self.Empuxo()))*1000000
       return s
    
    def N_p(self):
       n_p = (2*self.V_0()*self.Empuxo())/((self.A_0()**2)*(((1+self.F_o())*((self.V_9Ia_0())**(2)))-(m_0**2)))
       return n_p
    
    def N_t(self):
       n_ter = ((self.A_0()**2)*(((1+self.F_o())*((self.V_9Ia_0())**(2)))-(m_0**2)))/(2*self.F_o()*h_pr)
       return n_ter
    
    
    def N_o(self):
       n_o = self.N_t()*self.N_p()
       return n_o
    
def Nivel_motor(nivel):
    global t_t4, pi_dmax, pi_b, pi_n, e_c, e_t, n_b, pi_ab,n_ab, t_t7, e_f 
    if nivel == 0:           
        pi_dmax = 0.98
        e_c = 0.90
        e_f = 0.89  
        pi_b = 0.96
        n_b = 0.99
        e_t = 0.91
        pi_ab = 0.94
        n_ab = 0.95
        pi_n = 0.98   
        t_t4 = 3000 
        t_t7 = 3600
        if convert:
            t_t4 = t_t4/1.8 
            t_t7 = t_t7/1.8   

    elif nivel == 1:        
        pi_dmax = 0.88  #ajustavel
        e_t = 0.80  #ajustavel
        pi_n = 0.93  #ajustavel

        t_t4 = 1110 
        pi_b = 0.90
        e_c = 0.80
        n_b = 0.85

        e_f = 0.89
        t_t7 = 1390 
        pi_ab = 0.90
        n_ab = 0.85 
        if convert:
            t_t4 = 2000 
            t_t7 = 2500 
               

    elif nivel == 2:
        pi_dmax = 0.90  #ajustavel
        e_t = 0.83  #ajustavel
        pi_n = 0.93  #ajustavel

        t_t4 = 1390       
        pi_b = 0.92
        e_c = 0.84
        n_b = 0.91

        e_f = 0.89
        t_t7 = 1670
        pi_ab = 0.92 
        n_ab = 0.91
        if convert:
            t_t4 = 2500 
            t_t7 = 3000    

    elif nivel == 3:
        pi_dmax = 0.94  #ajustavel
        e_t = 0.87  #ajustavel
        pi_n = 0.95  #ajustavel

        t_t4 = 1780         
        pi_b = 0.94
        e_c = 0.88
        n_b = 0.98

        e_f = 0.89
        t_t7 = 2000 
        pi_ab = 0.94
        n_ab = 0.96
        if convert:
            t_t4 = 3200 
            t_t7 = 3600        
    
    elif nivel == 4:
        pi_dmax = 0.96  #ajustavel
        e_t = 0.89  #ajustavel
        pi_n = 0.97  #ajustavel

        t_t4 = 2000         
        pi_b = 0.95
        e_c = 0.90
        n_b = 0.99
        
        e_f = 0.89
        t_t7 = 2220
        pi_ab = 0.95
        n_ab = 0.99
        if convert:
            t_t4 = 3600
            t_t7 = 4000    

#####################################################################################

motor = Turbofan_mist_real()
nivel_motor = [3]

# Variaveis de entrada:
lista_pi_f = [0.63,0.99,1.5,2,7,12]
lista_pi_c = [0.01,40]

# Constantes de entrada em SI:
t_0 = 216.7         #[K]
c_pc = 1004     #[J/kg K]
c_pt = 1239    #[J/kg K]
c_pab = 1239    #[J/kg K]
h_pr = 42800000     #[J/kg]

gamma_c = 1.4 
gamma_t = 1.3
gamma_ab = 1.3 
n_m = 0.99
p_0Ip_9 = 1
m_6 = 0.5
m_0 = 1.8
pi_mmax = 0.95
pi_c_ = 17.4


after_burn = True

convert = True
if convert:
    # Constantes de entrada FORA do SI:
    t_0 = 390         #[R]
    c_pc = 0.240     #[Btu/lbm R]
    c_pt = 0.296    #[Btu/lbm R]
    c_pab = 0.296    #[Btu/lbm R]
    h_pr = 18400     #[Btu/lbm]
    
    #conversão
    t_0 = t_0/1.8             #[K]
    c_pc = (c_pc  * 4.1868)   #[kJ/kg K]
    c_pc  = c_pc * 1000        #[J/kg K]
    c_pt = (c_pt  * 4.1868)   #[kJ/kg K]
    c_pt = c_pt * 1000        #[J/kg K]
    c_pab = (c_pab  * 4.1868)   #[kJ/kg K]
    c_pab = c_pab * 1000        #[J/kg K]
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
    dado = format(t_t7, ".3f")
    print("%s\t%s\t%s" % ("| T_{t7}",dado,"K"))
    print("|-------------------------------------------|")
    dado = format(c_pc, ".3f")
    print("%s\t%s\t%s" % ("| c_{pc}",dado,"[J/kg K]"))
    print("|-------------------------------------------|")
    dado = format(c_pt, ".3f")
    print("%s\t%s\t%s" % ("| c_{pt}",dado,"[J/kg K]"))
    print("|-------------------------------------------|")
    dado = format(c_pab, ".3f")
    print("%s\t%s\t%s" % ("| c_{pab}",dado,"[J/kg K]"))
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
    dado = format(gamma_ab, ".3f")
    print("%s\t%s\t\t%s" % ("| gamma_{ab}",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_dmax, ".3f")
    print("%s\t%s\t\t%s" % ("| pi_{dmax}",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_b, ".3f")
    print("%s\t\t%s\t\t%s" % ("| pi_b",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_ab, ".3f")
    print("%s\t%s\t\t%s" % ("| pi_{ab}",dado," "))
    print("|-------------------------------------------|")
    dado = format(pi_mmax, ".3f")
    print("%s\t%s\t\t%s" % ("| pi_{Mmax}",dado," "))
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
    dado = format(n_ab, ".3f")
    print("%s\t%s\t\t%s" % ("| eta_{ab}",dado," "))
    print("|-------------------------------------------|")
    dado = format(n_m, ".3f")
    print("%s\t\t%s\t\t%s" % ("| eta_m",dado," "))
    print("|-------------------------------------------|")
    dado = format(p_0Ip_9, ".3f")
    print("%s\t%s\t\t%s" % ("| p_0/p_9",dado," "))
    print("|-------------------------------------------|")
    dado = format(m_6, ".3f")
    print("%s\t\t%s\t\t%s" % ("| m_6",dado," "))
    print("|-------------------------------------------|")
    print(" ")
    print("ººººººººººººººººººººººººººººººººººººººººººººº")

# Graficos:
plt.figure(figsize = ((15, 5)))
pi_c = np.linspace(min(lista_pi_c) ,max(lista_pi_c))  #Limites do eixo X

# Grafico 1 ----------------------------------------------------------
plt.subplot(2, 2, 1)

for j in range(len(lista_pi_f)):
    pi_f = lista_pi_f[j]
    for i in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.Empuxo(), label=r'$\pi_f=${},lvl={}'.format(pi_f,nivel_motor[i]))

y_ = 105*9.8067 
plt.axhline(y=y_, color='red', linestyle='--')
# x_ = 19.5
# plt.axvline(x=x_, color='red', linestyle='--')
# plt.scatter(x_, y_, color='red')
# plt.text(x_, y_, f'({x_:.1f}, {y_:.1f})', fontsize=12, ha='left', va='bottom')

            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(800, 1700)
plt.xlabel(r'$\pi_c$' , fontsize=15)
plt.ylabel(r'$F/ \dot m$   $[N / (kg /s)]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 2 ----------------------------------------------------------
plt.subplot(2, 2, 2)

for j in range(len(lista_pi_f)):
    pi_f = lista_pi_f[j]
    for i in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.S(), label=r'$\pi_f=${},lvl={}'.format(pi_f,nivel_motor[i]))

y_ = 1.845*28.325
plt.axhline(y=y_, color='red', linestyle='--')
# x_ = 17.3
# plt.axvline(x=x_, color='red', linestyle='--')
# plt.scatter(x_, y_, color='red')
# plt.text(x_, y_, f'({x_:.1f}, {y_:.1f})', fontsize=12, ha='left', va='bottom')
            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(40, 120)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$S$   $[(mg/s) / N]$ ' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 3 ----------------------------------------------------------
plt.subplot(2, 2, 3)

for j in range(len(lista_pi_f)):
    pi_f = lista_pi_f[j]
    for j in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.F(), label=r'$\pi_f=${},lvl={}'.format(pi_f,nivel_motor[i]))
            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$f$' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 4 ----------------------------------------------------------
plt.subplot(2, 2, 4)
for i in range(len(lista_pi_f)):
    pi_f = lista_pi_f[i]
    for i in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[i])
        plt.plot(pi_c, motor.Alpha(), label=r'$\pi_f=${},lvl={}'.format(pi_f,nivel_motor[i]))

y_ = 0.5
plt.axhline(y=y_, color='red', linestyle='--')
# x_ = 17.4
# plt.axvline(x=x_, color='red', linestyle='--')
# plt.scatter(x_, y_, color='red')
# plt.text(x_, y_, f'({x_:.1f}, {y_:.1f})', fontsize=12, ha='left', va='bottom')
            
plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
plt.ylim(0, 1.5)
plt.xlabel(r'$\pi_c$ ' , fontsize=15)
plt.ylabel(r'$\alpha *$' , fontsize=15)
plt.minorticks_on() # aparece a divisão
plt.grid()


# Plot Grafico ----------------------------------------------------------
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
plt.show()


    
# Tabela de saida:
print(" ")
for k in range(len(lista_pi_f)):
    pi_f = lista_pi_f[k]
    pi_c = pi_c_
    for j in range(len(nivel_motor)):
        Nivel_motor(nivel_motor[j])
        print(" \t\tTEC LEVEL = {}".format(nivel_motor[j]))
        print("|-------------------------------------------|")
        print("VALORES INICIAIS")
        dado = format(m_0, ".3f")
        print("%s\t\t%s\t\t%s" % ("| m_0",dado," "))
        print("|-------------------------------------------|")
        dado = format(pi_f, ".3f")
        print("%s\t\t%s\t\t%s" % ("| pi_f",dado," "))
        print("|-------------------------------------------|")
        print("%s\t\t%s\t\t%s" % ("| pi_c",pi_c," "))
        print("|-------------------------------------------|")
        print(" ")
        print(" ")
        print(" SAIDA\t\tVALOR\t\tUNIDADE")
        print("|-------------------------------------------|")
        dado = format(motor.S(), ".3f")
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
        dado = format(motor.F_AB(), ".3f")
        print("%s\t%s\t%s" % ("| f_{AB}",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.F_o(), ".3f")
        print("%s\t\t%s\t%s" % ("| f_O",dado," "))
        print("|-------------------------------------------|")
        dado = format(motor.Alpha(), ".3f")
        print("%s\t\t%s\t%s" % ("| alpha",dado," "))
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
        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        print(" ")




