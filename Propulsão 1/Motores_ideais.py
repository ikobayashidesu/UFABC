import matplotlib.pyplot as plt
import numpy as np

turbofan_mix_alpha_ = False
turbofan_mix_pi_f_ = False
turbofan_mix_alpha_AF_ = False
turbofan_mix_pi_f_AF_ = True

class Turbofan_mix_alpha:
    def __init__(self):
        s = 1

    def R_(self):
        R_ = ((gamma-1)/(gamma))*c_p
        self.R = R_
        return R_

    def A_0(self):
        a_0_ = (gamma*self.R_()*t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_
    
    def Tau_de_r(self):
        tau_r_ = 1 + (((gamma - 1)/(2))*(m_0**2))
        self.tau_r = tau_r_
        return tau_r_
    
    def Tau_de_lambida(self):
        tau_lambida_ = t_t4/t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_
    
    def Tau_de_c(self):
        tau_c_ = (pi_c)**((gamma-1)/(gamma))
        self.tau_c = tau_c_
        return tau_c_
    
    def Tau_de_f(self):
        tau_f_ = (((self.Tau_de_lambida())/(self.Tau_de_r()))-(self.Tau_de_c()-1)+ alpha)/(((self.Tau_de_lambida())/(self.Tau_de_r()*self.Tau_de_c()))+ alpha)
        self.tau_f = tau_f_
        return tau_f_
    
    def Pi_de_f(self):
        pi_f_ = (self.Tau_de_f())**((gamma)/(gamma-1))
        self.pi_f = pi_f_
        return pi_f_
    
    def Tau_de_t(self):
        tau_t_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1+(alpha*(self.Tau_de_f()-1))))
        self.tau_t = tau_t_
        return tau_t_
    
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Tau_de_m(self):
        tau_m_ = ((1)/(1+alpha))*(1+(alpha*((self.Tau_de_r()*self.Tau_de_f())/(self.Tau_de_lambida()*self.Tau_de_t()))))
        self.tau_m = tau_m_
        return tau_m_
    
    def Razao_de_temperatura(self):
        t9t0_ = (self.Tau_de_lambida()*self.Tau_de_t()*self.Tau_de_m())/(self.Tau_de_r()*self.Tau_de_f()) 
        self.t9t0 = t9t0_
        return t9t0_
    
    def Mach_9(self):
        m9_ = (((2)/(gamma - 1))*((self.Tau_de_r()*self.Tau_de_f())-1))**(1/2)
        self.m9 = m9_
        return m9_

    def Razao_de_velocidade(self):
        raz_vel_ = ((self.Razao_de_temperatura())**(1/2))*(self.Mach_9())
        self.raz_vel = raz_vel_ 
        return raz_vel_ 
    
    def Razao_de_combustivel_ar_o(self):
        fo_ = (self.Razao_combustivel_ar())/(1+alpha)
        self.fo = fo_
        return fo_

    def Empuxo(self):
        empuxo_ = self.A_0()*(self.Razao_de_velocidade()-m_0)
        self.empuxo = empuxo_
        return empuxo_
   
    def Consumo_especifico(self):
        s_ = self.Razao_de_combustivel_ar_o()/self.Empuxo()
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = ((gamma-1)/(2))*((c_p*t_0)/(self.Razao_de_combustivel_ar_o()*h_pr))*((self.Razao_de_velocidade()**2)-(m_0**2))
        self.n_t = n_t_
        return n_t_
    
    def Eficiencia_propulsiva(self):
        n_p_ = (2*m_0)/(self.Razao_de_velocidade() + m_0)
        self.n_p = n_p_
        return n_p_
    
    def Eficiencia_total(self):
        n_o_ = self.Eficiencia_termica()*self.Eficiencia_propulsiva()
        self.n_o = n_o_ 
        return n_o_

class Turbofan_mix_pi_f:
    def __init__(self):
        s = 1

    def R_(self):
        R_ = ((gamma-1)/(gamma))*c_p
        self.R = R_
        return R_

    def A_0(self):
        a_0_ = (gamma*self.R_()*t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_
    
    def Tau_de_r(self):
        tau_r_ = 1 + (((gamma - 1)/(2))*(m_0**2))
        self.tau_r = tau_r_
        return tau_r_
    
    def Tau_de_lambida(self):
        tau_lambida_ = t_t4/t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_
    
    def Tau_de_c(self):
        tau_c_ = (pi_c)**((gamma-1)/(gamma))
        self.tau_c = tau_c_
        return tau_c_
    
    def Tau_de_f(self):
        tau_f_ = pi_f**((gamma-1)/(gamma))
        self.tau_f = tau_f_
        return tau_f_
    
    def Alpha(self):
        alpha_ = ((self.Tau_de_lambida()*(self.Tau_de_c()-self.Tau_de_f()))/(self.Tau_de_r()*self.Tau_de_c()*(self.Tau_de_f()-1)))-((self.Tau_de_c()-1)/(self.Tau_de_f()-1))
        self.alpha = alpha_
        return alpha_
    
    def Tau_de_t(self):
        tau_t_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1+(self.Alpha()*(self.Tau_de_f()-1))))
        self.tau_t = tau_t_
        return tau_t_
    
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Tau_de_m(self):
        tau_m_ = ((1)/(1+self.Alpha()))*(1+(self.Alpha()*((self.Tau_de_r()*self.Tau_de_f())/(self.Tau_de_lambida()*self.Tau_de_t()))))
        self.tau_m = tau_m_
        return tau_m_
    
    def Razao_de_temperatura(self):
        t9t0_ = (self.Tau_de_lambida()*self.Tau_de_t()*self.Tau_de_m())/(self.Tau_de_r()*self.Tau_de_f()) 
        self.t9t0 = t9t0_
        return t9t0_
    
    def Mach_9(self):
        m9_ = (((2)/(gamma - 1))*((self.Tau_de_r()*self.Tau_de_f())-1))**(1/2)
        self.m9 = m9_
        return m9_

    def Razao_de_velocidade(self):
        raz_vel_ = ((self.Razao_de_temperatura())**(1/2))*(self.Mach_9())
        self.raz_vel = raz_vel_ 
        return raz_vel_ 
    
    def Razao_de_combustivel_ar_o(self):
        fo_ = (self.Razao_combustivel_ar())/(1+self.Alpha())
        self.fo = fo_
        return fo_

    def Empuxo(self):
        empuxo_ = self.A_0()*(self.Razao_de_velocidade()-m_0)
        self.empuxo = empuxo_
        return empuxo_
   
    def Consumo_especifico(self):
        s_ = self.Razao_de_combustivel_ar_o()/self.Empuxo()
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = ((gamma-1)/(2))*((c_p*t_0)/(self.Razao_de_combustivel_ar_o()*h_pr))*((self.Razao_de_velocidade()**2)-(m_0**2))
        self.n_t = n_t_
        return n_t_
    
    def Eficiencia_propulsiva(self):
        n_p_ = (2*m_0)/(self.Razao_de_velocidade() + m_0)
        self.n_p = n_p_
        return n_p_
    
    def Eficiencia_total(self):
        n_o_ = self.Eficiencia_termica()*self.Eficiencia_propulsiva()
        self.n_o = n_o_ 
        return n_o_

class Turbofan_mix_alpha_AF:
    def __init__(self):
        s = 1

    def R_(self):
        R_ = ((gamma-1)/(gamma))*c_p
        self.R = R_
        return R_

    def A_0(self):
        a_0_ = (gamma*self.R_()*t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_
    
    def Tau_de_r(self):
        tau_r_ = 1 + (((gamma - 1)/(2))*(m_0**2))
        self.tau_r = tau_r_
        return tau_r_
    
    def Tau_de_lambida(self):
        tau_lambida_ = t_t4/t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_
    
    def Tau_de_c(self):
        tau_c_ = (pi_c)**((gamma-1)/(gamma))
        self.tau_c = tau_c_
        return tau_c_
    
    def Tau_de_f(self):
        tau_f_ = (((self.Tau_de_lambida())/(self.Tau_de_r()))-(self.Tau_de_c()-1)+ alpha)/(((self.Tau_de_lambida())/(self.Tau_de_r()*self.Tau_de_c()))+ alpha)
        self.tau_f = tau_f_
        return tau_f_
    
    def Pi_de_f(self):
        pi_f_ = (self.Tau_de_f())**((gamma)/(gamma-1))
        self.pi_f = pi_f_
        return pi_f_
    
    def Tau_de_t(self):
        tau_t_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1+(alpha*(self.Tau_de_f()-1))))
        self.tau_t = tau_t_
        return tau_t_
    
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Tau_de_m(self):
        tau_m_ = ((1)/(1+alpha))*(1+(alpha*((self.Tau_de_r()*self.Tau_de_f())/(self.Tau_de_lambida()*self.Tau_de_t()))))
        self.tau_m = tau_m_
        return tau_m_
    
    def Tau_ab(self):
        tau_ab_ = t_t7/t_0
        self.tau_ab = tau_ab_
        return tau_ab_
    
    def Razao_de_combustivel_ar_AB(self):
        f_ab_ = ((c_p*t_0)/(h_pr))*(self.Tau_ab()-(self.Tau_de_lambida()*self.Tau_de_t()*self.Tau_de_m()))
        self.f_ab = f_ab_
        return f_ab_
    
    def Razao_de_temperatura(self):
        t9t0_ = (self.Tau_ab())/(self.Tau_de_r()*self.Tau_de_f())
        self.t9t0 = t9t0_
        return t9t0_
    
    def Mach_9(self):
        m9_ = (((2)/(gamma - 1))*((self.Tau_de_r()*self.Tau_de_f())-1))**(1/2)
        self.m9 = m9_
        return m9_

    def Razao_de_velocidade(self):
        raz_vel_ = ((self.Razao_de_temperatura())**(1/2))*(self.Mach_9())
        self.raz_vel = raz_vel_ 
        return raz_vel_ 
    
    def Razao_de_combustivel_ar_o(self):
        fo_ = ((self.Razao_combustivel_ar())/(1+alpha)) + self.Razao_de_combustivel_ar_AB()
        self.fo = fo_
        return fo_

    def Empuxo(self):
        empuxo_ = self.A_0()*(self.Razao_de_velocidade()-m_0)
        self.empuxo = empuxo_
        return empuxo_
   
    def Consumo_especifico(self):
        s_ = self.Razao_de_combustivel_ar_o()/self.Empuxo()
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = ((gamma-1)/(2))*((c_p*t_0)/(self.Razao_de_combustivel_ar_o()*h_pr))*((self.Razao_de_velocidade()**2)-(m_0**2))
        self.n_t = n_t_
        return n_t_
    
    def Eficiencia_propulsiva(self):
        n_p_ = (2*m_0)/(self.Razao_de_velocidade() + m_0)
        self.n_p = n_p_
        return n_p_
    
    def Eficiencia_total(self):
        n_o_ = self.Eficiencia_termica()*self.Eficiencia_propulsiva()
        self.n_o = n_o_ 
        return n_o_

class Turbofan_mix_pi_f_AF:
    def __init__(self):
        s = 1

    def R_(self):
        R_ = ((gamma-1)/(gamma))*c_p
        self.R = R_
        return R_

    def A_0(self):
        a_0_ = (gamma*self.R_()*t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_
    
    def Tau_de_r(self):
        tau_r_ = 1 + (((gamma - 1)/(2))*(m_0**2))
        self.tau_r = tau_r_
        return tau_r_
    
    def Tau_de_lambida(self):
        tau_lambida_ = t_t4/t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_
    
    def Tau_de_c(self):
        tau_c_ = (pi_c)**((gamma-1)/(gamma))
        self.tau_c = tau_c_
        return tau_c_
    
    def Tau_de_f(self):
        tau_f_ = pi_f**((gamma-1)/(gamma))
        self.tau_f = tau_f_
        return tau_f_
    
    def Alpha(self):
        alpha_ = ((self.Tau_de_lambida()*(self.Tau_de_c()-self.Tau_de_f()))/(self.Tau_de_r()*self.Tau_de_c()*(self.Tau_de_f()-1)))-((self.Tau_de_c()-1)/(self.Tau_de_f()-1))
        self.alpha = alpha_
        return alpha_
    
    def Tau_de_t(self):
        tau_t_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1+(self.Alpha()*(self.Tau_de_f()-1))))
        self.tau_t = tau_t_
        return tau_t_
    
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Tau_de_m(self):
        tau_m_ = ((1)/(1+self.Alpha()))*(1+(self.Alpha()*((self.Tau_de_r()*self.Tau_de_f())/(self.Tau_de_lambida()*self.Tau_de_t()))))
        self.tau_m = tau_m_
        return tau_m_
    
    def Tau_ab(self):
        tau_ab_ = t_t7/t_0
        self.tau_ab = tau_ab_
        return tau_ab_
    
    def Razao_de_combustivel_ar_AB(self):
        f_ab_ = ((c_p*t_0)/(h_pr))*(self.Tau_ab()-(self.Tau_de_lambida()*self.Tau_de_t()*self.Tau_de_m()))
        self.f_ab = f_ab_
        return f_ab_
    
    def Razao_de_temperatura(self):
        t9t0_ = (self.Tau_ab())/(self.Tau_de_r()*self.Tau_de_f())
        self.t9t0 = t9t0_
        return t9t0_
    
    def Mach_9(self):
        m9_ = (((2)/(gamma - 1))*((self.Tau_de_r()*self.Tau_de_f())-1))**(1/2)
        self.m9 = m9_
        return m9_

    def Razao_de_velocidade(self):
        raz_vel_ = ((self.Razao_de_temperatura())**(1/2))*(self.Mach_9())
        self.raz_vel = raz_vel_ 
        return raz_vel_ 

    def Razao_de_combustivel_ar_o(self):
        fo_ = ((self.Razao_combustivel_ar())/(1+self.Alpha()))+self.Razao_de_combustivel_ar_AB()
        self.fo = fo_
        return fo_

    def Empuxo(self):
        empuxo_ = self.A_0()*(self.Razao_de_velocidade()-m_0)
        self.empuxo = empuxo_
        return empuxo_
    
    def Consumo_especifico(self):
        s_ = self.Razao_de_combustivel_ar_o()/self.Empuxo()
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = ((gamma-1)/(2))*((c_p*t_0)/(self.Razao_de_combustivel_ar_o()*h_pr))*((self.Razao_de_velocidade()**2)-(m_0**2))
        self.n_t = n_t_
        return n_t_
    
    def Eficiencia_propulsiva(self):
        n_p_ = (2*m_0)/(self.Razao_de_velocidade() + m_0)
        self.n_p = n_p_
        return n_p_
    
    def Eficiencia_total(self):
        n_o_ = self.Eficiencia_termica()*self.Eficiencia_propulsiva()
        self.n_o = n_o_ 
        return n_o_
###################################################################################################################################################################

if turbofan_mix_alpha_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # SEM PÓS COMBUSTOR
    # COM RAZÃO BYPESS

    motor_2 = Turbofan_mix_alpha()

    # Variaveis de entrada:
    m_0 = 2.5
    lista_alpha = [0.5,1]
    pi_c = 16
    lista_t_t4 = [1600,2200] #[K]

    # Constantes de entrada:
    gamma = 1.4 
    t_0 = 216.7         #[K]      
    c_p = 1.004    #[kJ/kg K]
    h_pr = 42800   #[kJ/kg]

    #conversões:
    c_p  = c_p * 1000              #[J/kg K]
    h_pr = h_pr * 1000              #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Pi_de_f(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$\pi_f$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Empuxo(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Consumo_especifico()*1000000, label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$S [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_alpha)):
            alpha = lista_alpha[j]
            print(f"Razão de Bypass (alpha)________________ {alpha} ")
            dado = format(motor_2.Pi_de_f(), ".3f")
            print(f"Razão de pressão do fan (PI_f)_________ {dado} ")
            dado = format(motor_2.Empuxo(), ".3f")
            print(f"Empuxo (F)_____________________________ {dado} [N/(kg s)]")
            dado = format(motor_2.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)___ {dado} [(mg s)/(N)]")
            print(" ")
        print("**************************************************")

if turbofan_mix_pi_f_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # SEM PÓS COMBUSTOR
    # COM RAZÃO DE PRESSÃO DO FAN

    motor = Turbofan_mix_pi_f()

    # Variaveis de entrada:
    m_0 = 2
    lista_pi_f = [2,5]
    pi_c = 20
    lista_t_t4 = [3000,4000] #[R]

    # Constantes de entrada:
    gamma = 1.4 
    t_0 = 390         #[R]      
    c_p = 0.24    #[Btu/lbm R]
    h_pr = 18400   #[Btu/lbm]

    #conversões:
    t_0 = t_0/1.8                 #[K]
    for i in range(len(lista_t_t4)):
      lista_t_t4[i] = lista_t_t4[i]/1.8       #[K]
    c_p  = (c_p  * 4.1868)      #[kJ/kg K]
    c_p  = c_p * 1000              #[J/kg K]
    h_pr = h_pr * 2.326           #[kJ/kg]
    h_pr = h_pr * 1000              #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Alpha(), label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$\alpha$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Empuxo(), label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Consumo_especifico()*1000000, label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$S [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_pi_f)):
            pi_f = lista_pi_f[j]
            print(f"Razão de pressão do fan (PI_f)_________ {pi_f} ")
            dado = format(motor.Alpha(), ".3f")
            print(f"Razão Bypass (alpha)___________________ {dado} ")
            dado = format(motor.Empuxo(), ".3f")
            print(f"Empuxo (F)_____________________________ {dado} [N/(kg s)]")
            dado = format(motor.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)___ {dado} [(mg s)/(N)]")
            print(" ")
        print("**************************************************")

if turbofan_mix_alpha_AF_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # COM PÓS COMBUSTOR
    # COM RAZÃO BYPESS

    motor_3 = Turbofan_mix_alpha_AF()

    # Variaveis de entrada:
    m_0 = 2.5
    lista_alpha = [0.5,1]
    pi_c = 16
    lista_t_t4 = [1600,2200] #[K]

    # Constantes de entrada:
    gamma = 1.4 
    t_0 = 216.7         #[K]
    t_t7 = 2200         #[K]        
    c_p = 1.004    #[kJ/kg K]
    h_pr = 42800   #[kJ/kg]

    #conversões:
    c_p  = c_p * 1000              #[J/kg K]
    h_pr = h_pr * 1000              #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Pi_de_f(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$\pi_f$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Empuxo(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Consumo_especifico()*1000000, label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$S [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_alpha)):
            alpha = lista_alpha[j]
            print(f"Razão de Bypass (alpha)________________ {alpha} ")
            dado = format(motor_3.Pi_de_f(), ".3f")
            print(f"Razão de pressão do fan (PI_f)_________ {dado} ")
            dado = format(motor_3.Empuxo(), ".3f")
            print(f"Empuxo (F)_____________________________ {dado} [N/(kg s)]")
            dado = format(motor_3.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)___ {dado} [(mg s)/(N)]")
            print(" ")
        print("**************************************************")
    
if turbofan_mix_pi_f_AF_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # COM PÓS COMBUSTOR
    # COM RAZÃO DE PRESSÃO DO FAN

    motor_4 = Turbofan_mix_pi_f_AF()

    # Variaveis de entrada:
    m_0 = 2
    lista_pi_f = [2,5]
    pi_c = 20
    lista_t_t4 = [3000, 4000] #[R]

    # Constantes de entrada:
    gamma = 1.4
    t_t7 = 4000       #[R]
    t_0 = 390         #[R]      
    c_p = 0.24    #[Btu/lbm R]
    h_pr = 18400   #[Btu/lbm]

    #conversões:
    t_t7 = t_t7/1.8              #[K]
    t_0 = t_0/1.8                 #[K]
    for i in range(len(lista_t_t4)):
      lista_t_t4[i] = lista_t_t4[i]/1.8       #[K]
    c_p  = (c_p  * 4.1868)      #[kJ/kg K]
    c_p  = c_p * 1000              #[J/kg K]
    h_pr = h_pr * 2.326           #[kJ/kg]
    h_pr = h_pr * 1000              #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor_4.Alpha(), label=r'$\pi_f = ${}'.format(pi_f))
    
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$\alpha$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor_4.Consumo_especifico(), label=r'$\pi_f = ${}'.format(pi_f)) 
    
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor_4.Consumo_especifico()*1000000, label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4} [K]$ ' , fontsize=15)
    plt.ylabel(r'$S [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_pi_f)):
            pi_f = lista_pi_f[j]
            print(f"Razão de pressão do fan (PI_f)_________ {pi_f} ")
            dado = format(motor_4.Alpha(), ".3f")
            print(f"Razão Bypass (alpha)___________________ {dado} ")
            dado = format(motor_4.Empuxo(), ".3f")
            print(f"Empuxo (F)_____________________________ {dado} [N/(kg s)]")
            dado = format(motor_4.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)___ {dado} [(mg s)/(N)]")
            print(" ")
        print("**************************************************")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------


