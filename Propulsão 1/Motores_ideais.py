import matplotlib.pyplot as plt
import numpy as np

turbofan_mix_alpha_ = False
turbofan_mix_pi_f_ = False
turbofan_mix_alpha_AF_ = False
turbofan_mix_pi_f_AF_ = False
turboprop = False

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

class Turboprop:
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
     
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Tau_de_th(self):
        tau_th_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1))
        self.tau_th = tau_th_
        return tau_th_
    
    def Tau_de_tl(self):
        tau_tl_= (t_t)/(self.Tau_de_th())
        self.tau_de_tl = tau_tl_
        return tau_tl_
    
    def Razao_de_velocidade(self):
        raz_vel_ = (((2)/(gamma-1))*((self.Tau_de_lambida()*t_t)-((self.Tau_de_lambida())/(self.Tau_de_r()*self.Tau_de_c()))))**(1/2)
        self.raz_vel = raz_vel_ 
        return raz_vel_ 
    
    def C_c(self):
        c_c_ = (gamma-1)*(m_0)*(self.Razao_de_velocidade()-m_0)
        self.c_c = c_c_
        return c_c_
    
    def C_de_prop(self):
        c_prop_ = n_prop*self.Tau_de_lambida()*self.Tau_de_th()*(1-self.Tau_de_tl())
        self.c_prop = c_prop_
        return c_prop_
    
    def C_de_total(self):
        c_total_ = self.C_de_prop() + self.C_c()
        self.c_total = c_total_
        return c_total_

    def Empuxo(self):
        empuxo_ = (self.C_de_total()*c_p*t_0)/(self.A_0()*m_0)
        self.empuxo = empuxo_
        return empuxo_
   
    def Consumo_especifico(self):
        s_ = (self.Razao_combustivel_ar()/self.Empuxo())*1000000
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = 1 - ((1)/(self.Tau_de_r()*self.Tau_de_c()))
        self.n_t = n_t_
        return n_t_
    
    def Eficiencia_propulsiva(self):
        n_p_ = self.Eficiencia_total()/self.Eficiencia_termica()
        self.n_p = n_p_
        return n_p_
    
    def Eficiencia_total(self):
        n_o_ = (self.C_de_total())/(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.n_o = n_o_ 
        return n_o_
    
    def Tau_t_otimo(self):
        tau_t_otimo_ = ((1)/(self.Tau_de_r()*self.Tau_de_c()))+((((gamma-1)/(2))*(m_0**2))/(self.Tau_de_lambida()*(n_prop**2)))
        self.tau_t_otimo = tau_t_otimo_
        return tau_t_otimo_
###################################################################################################################################################################

if turbofan_mix_alpha_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # SEM PÓS COMBUSTOR
    # COM RAZÃO BYPESS

    motor_2 = Turbofan_mix_alpha()

    # Variaveis de entrada:
    lista_alpha = [0.5,1]
    lista_t_t4 = [1600,2200]  #[K]

    # Constantes de entrada:
    pi_c = 16
    m_0 = 2.5
    gamma = 1.4 
    t_0 = 216.7    #[K]      
    c_p = 1.004    #[kJ/kg K]
    h_pr = 42800   #[kJ/kg]

    #conversões:
    c_p  = c_p * 1000    #[J/kg K]
    h_pr = h_pr * 1000   #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))
    t_t4_ = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Pi_de_f(), label=r'$\alpha = ${}'.format(alpha))
           
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$\pi_f$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Empuxo(), label=r'$\alpha = ${}'.format(alpha))
             
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $[K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m$   $ [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_2.Consumo_especifico()*1000000, label=r'$\alpha = ${}'.format(alpha))
               
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$S$   $ [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 4 ----------------------------------------------------------
    plt.subplot(2, 2, 4)
    lista_S = []
    lista_F = []
    for j in range(len(lista_alpha)):
        alpha = lista_alpha[j]
        for i in range(len(t_t4_)):
            t_t4 = t_t4_[i]
            lista_S.append(motor_2.Consumo_especifico()*1000000)
            lista_F.append(motor_2.Empuxo())
        plt.plot(lista_F, lista_S, label=r'$\alpha = ${}'.format(alpha))
        lista_F.clear()
        lista_S.clear()

         
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$F/ \dot m$   $ [N / (kg /s)]$' , fontsize=15)
    plt.ylabel(r'$S$   $ [(mg/s) / N]$ ' , fontsize=15)
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

            dado = format(motor_2.Empuxo(), ".3f")
            print(f"Empuxo (F)_____________________________ {dado} [N/(kg s)]")

            dado = format(motor_2.Razao_combustivel_ar(), ".3f")
            print(f"Razão combustivel ar (f)_______________ {dado}")

            dado = format(motor_2.Razao_de_combustivel_ar_o(), ".3f")
            print(f"Razão combustivel ar total (f_0)_______ {dado}")

            dado = format(motor_2.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)___ {dado} [(mg s)/(N)]")

            dado = format(motor_2.Eficiencia_termica(), ".3f")
            print(f"Eficiencia termica (n_t)________________ {dado} ")

            dado = format(motor_2.Eficiencia_propulsiva(), ".3f")
            print(f"Eficiencia propulsiva (n_p)_____________ {dado} ")

            dado = format(motor_2.Eficiencia_termica(), ".3f")
            print(f"Eficiencia total (n_o)__________________ {dado} ")

            dado = format(motor_2.Pi_de_f(), ".3f")
            print(f"Razão de pressão do fan (PI_f)_________ {dado} ")

            print(" ")
        print("**************************************************")

if turbofan_mix_pi_f_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # SEM PÓS COMBUSTOR
    # COM RAZÃO DE PRESSÃO DO FAN

    motor = Turbofan_mix_pi_f()

    # Variaveis de entrada:
    lista_pi_f = [2,5]
    lista_t_t4 = [3000,4000] #[R]

    # Constantes de entrada:
    pi_c = 20
    m_0 = 2
    gamma = 1.4 
    t_0 = 390      #[R]      
    c_p = 0.24     #[Btu/lbm R]
    h_pr = 18400   #[Btu/lbm]

    #conversões:
    for i in range(len(lista_t_t4)):
      lista_t_t4[i] = lista_t_t4[i]/1.8  #[K]
    t_0 = t_0/1.8                        #[K]
    c_p  = (c_p  * 4.1868)               #[kJ/kg K]
    c_p  = c_p * 1000                    #[J/kg K]
    h_pr = h_pr * 2.326                  #[kJ/kg]
    h_pr = h_pr * 1000                   #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Alpha(), label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $[K]$ ' , fontsize=15)
    plt.ylabel(r'$\alpha$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Empuxo(), label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $[K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m$   $[N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor.Consumo_especifico()*1000000, label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $[K]$ ' , fontsize=15)
    plt.ylabel(r'$S$   $[(mg/s) / N]$ ' , fontsize=15)
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
            print(f"Razão de pressão do fan (PI_f)__________ {pi_f} ")
            
            dado = format(motor.Empuxo(), ".3f")
            print(f"Empuxo (F)______________________________ {dado} [N/(kg s)]")

            dado = format(motor.Razao_combustivel_ar(), ".3f")
            print(f"Razão combustivel ar (f)________________ {dado}")

            dado = format(motor.Razao_de_combustivel_ar_o(), ".3f")
            print(f"Razão combustivel ar total (f_0)________ {dado}")

            dado = format(motor.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)____ {dado} [(mg s)/(N)]")

            dado = format(motor.Eficiencia_termica(), ".3f")
            print(f"Eficiencia termica (n_t)________________ {dado} ")

            dado = format(motor.Eficiencia_propulsiva(), ".3f")
            print(f"Eficiencia propulsiva (n_p)_____________ {dado} ")

            dado = format(motor.Eficiencia_termica(), ".3f")
            print(f"Eficiencia total (n_o)__________________ {dado} ")

            dado = format(motor.Alpha(), ".3f")
            print(f"Razão Bypass (alpha)____________________ {dado} ")
            print(" ")
        print("**************************************************")

if turbofan_mix_alpha_AF_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # COM PÓS COMBUSTOR
    # COM RAZÃO BYPESS

    motor_3 = Turbofan_mix_alpha_AF()

    # Variaveis de entrada:
    lista_alpha = [0.5,1]
    lista_t_t4 = [1600,2200] #[K]

    # Constantes de entrada:
    pi_c = 16
    m_0 = 2.5
    gamma = 1.4 
    t_0 = 216.7    #[K]
    t_t7 = 2200    #[K]        
    c_p = 1.004    #[kJ/kg K]
    h_pr = 42800   #[kJ/kg]

    #conversões:
    c_p  = c_p * 1000    #[J/kg K]
    h_pr = h_pr * 1000   #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X
    t_t4_ = np.linspace(min(lista_t_t4) ,max(lista_t_t4)) 

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Pi_de_f(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$\pi_f$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Empuxo(), label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m$   $ [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    plt.subplot(2, 2, 3)
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        plt.plot(t_t4, motor_3.Consumo_especifico()*1000000, label=r'$\alpha = ${}'.format(alpha))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$S$   $ [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 4 ----------------------------------------------------------
    plt.subplot(2, 2, 4)
    lista_S = []
    lista_F = []
    lista_F.clear()
    lista_S.clear()
    for j in range(len(lista_alpha)):
        alpha = lista_alpha[j]
        for i in range(len(t_t4_)):
            t_t4 = t_t4_[i]
            lista_S.append(motor_3.Consumo_especifico()*1000000)
            lista_F.append(motor_3.Empuxo())
        plt.plot(lista_F, lista_S, label=r'$\alpha = ${}'.format(alpha))
        lista_F.clear()
        lista_S.clear()

         
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$F/ \dot m$   $ [N / (kg /s)]$' , fontsize=15)
    plt.ylabel(r'$S$   $ [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)___________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_alpha)):
            alpha = lista_alpha[j]
            print(f"Razão de Bypass (alpha)___________________ {alpha} ")
            
            dado = format(motor_3.Empuxo(), ".3f")
            print(f"Empuxo (F)________________________________ {dado} [N/(kg s)]")

            dado = format(motor_3.Razao_combustivel_ar(), ".3f")
            print(f"Razão combustivel ar (f)__________________ {dado}")

            dado = format(motor_3.Razao_de_combustivel_ar_AB(), ".3f")
            print(f"Razão combustivel ar Pós combustor (f_ab)_ {dado}")

            dado = format(motor_3.Razao_de_combustivel_ar_o(), ".3f")
            print(f"Razão combustivel ar total (f_0)__________ {dado}")

            dado = format(motor_3.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)______ {dado} [(mg s)/(N)]")

            dado = format(motor_3.Eficiencia_termica(), ".3f")
            print(f"Eficiencia termica (n_t)__________________ {dado} ")

            dado = format(motor_3.Eficiencia_propulsiva(), ".3f")
            print(f"Eficiencia propulsiva (n_p)_______________ {dado} ")

            dado = format(motor_3.Eficiencia_termica(), ".3f")
            print(f"Eficiencia total (n_o)____________________ {dado} ")

            dado = format(motor_3.Pi_de_f(), ".3f")
            print(f"Razão de pressão do fan (PI_f)____________ {dado} ")

            print(" ")
        print("**************************************************")
    
if turbofan_mix_pi_f_AF_:

    # TRBOFAN COM ESCOAMENTO MISTURADO 
    # COM PÓS COMBUSTOR
    # COM RAZÃO DE PRESSÃO DO FAN

    motor_4 = Turbofan_mix_pi_f_AF()

    # Variaveis de entrada:
    lista_pi_f = [2,5]
    lista_t_t4 = [3000, 4000] #[R]

    # Constantes de entrada:
    pi_c = 20
    m_0 = 2
    gamma = 1.4
    t_t7 = 4000       #[R]
    t_0 = 390         #[R]      
    c_p = 0.24        #[Btu/lbm R]
    h_pr = 18400      #[Btu/lbm]

    #conversões:
    for i in range(len(lista_t_t4)):
      lista_t_t4[i] = lista_t_t4[i]/1.8  #[K]
    t_t7 = t_t7/1.8                      #[K]
    t_0 = t_0/1.8                        #[K]
    c_p  = (c_p  * 4.1868)               #[kJ/kg K]
    c_p  = c_p * 1000                    #[J/kg K]
    h_pr = h_pr * 2.326                  #[kJ/kg]
    h_pr = h_pr * 1000                   #[J/kg]
    
    # Tabelas:
    plt.figure(figsize = ((12, 6)))
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X
    t_t4_ = np.linspace(min(lista_t_t4) ,max(lista_t_t4)) 

    # Grafico 1 ----------------------------------------------------------
    plt.subplot(2, 2, 1)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor_4.Alpha(), label=r'$\pi_f = ${}'.format(pi_f))
    
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$\alpha$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 2 ----------------------------------------------------------
    plt.subplot(2, 2, 2)
    
    # for i in range(len(lista_pi_f)):
    #     pi_f = lista_pi_f[i]
    #     plt.plot(t_t4, motor_4.Consumo_especifico(), label=r'$\pi_f = ${}'.format(pi_f)) 
    
    lista_F = []
    lista_t4 = []
    for j in range(len(lista_pi_f)):
        pi_f = lista_pi_f[j]
        for i in range(len(t_t4_)):
            t_t4 = t_t4_[i]
            lista_t4.append(t_t4_[i])
            lista_F.append(motor_4.Empuxo())
        plt.plot(lista_t4, lista_F, label=r'$\pi_f = ${}'.format(pi_f))
        lista_F.clear()
        lista_t4.clear()

    
    
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$F/ \dot m$   $ [N / (kg /s)]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Grafico 3 ----------------------------------------------------------
    t_t4 = np.linspace(min(lista_t_t4) ,max(lista_t_t4))  #Limites do eixo X
    plt.subplot(2, 2, 3)
    for i in range(len(lista_pi_f)):
        pi_f = lista_pi_f[i]
        plt.plot(t_t4, motor_4.Consumo_especifico()*1000000, label=r'$\pi_f = ${}'.format(pi_f))
            
    plt.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$T_{t4}$   $ [K]$ ' , fontsize=15)
    plt.ylabel(r'$S$   $ [(mg/s) / N]$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
    plt.show()
    
    for i in range(len(lista_t_t4)):
        t_t4 = lista_t_t4[i]
        print(" ")
        t_t4 = round(t_t4, 3)
        print(f"Temperatura no combustor (T_t4)___________ {t_t4 } [K]")
        print(" ")
        for j in range(len(lista_pi_f)):
            pi_f = lista_pi_f[j]
            print(f"Razão de pressão do fan (PI_f)____________ {pi_f} ")

            dado = format(motor_4.Empuxo(), ".3f")
            print(f"Empuxo (F)________________________________ {dado} [N/(kg s)]")

            dado = format(motor_4.Razao_combustivel_ar(), ".3f")
            print(f"Razão combustivel ar (f)__________________ {dado}")

            dado = format(motor_4.Razao_de_combustivel_ar_AB(), ".3f")
            print(f"Razão combustivel ar Pós combustor (f_ab)_ {dado}")

            dado = format(motor_4.Razao_de_combustivel_ar_o(), ".3f")
            print(f"Razão combustivel ar total (f_0)__________ {dado}")

            dado = format(motor_4.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)______ {dado} [(mg s)/(N)]")

            dado = format(motor_4.Eficiencia_termica(), ".3f")
            print(f"Eficiencia termica (n_t)__________________ {dado} ")

            dado = format(motor_4.Eficiencia_propulsiva(), ".3f")
            print(f"Eficiencia propulsiva (n_p)_______________ {dado} ")

            dado = format(motor_4.Eficiencia_termica(), ".3f")
            print(f"Eficiencia total (n_o)____________________ {dado} ")

            dado = format(motor_4.Alpha(), ".3f")
            print(f"Razão Bypass (alpha)______________________ {dado} ")
            print(" ")
        print("**************************************************")

if turboprop:
    motor_5 = Turboprop()

    # Variaveis de entrada:
    m_0 = 0.65
    lista_pi_c = [2,40]
    lista_t_t = [0.8,0.7,0.6,0.5,0.4]

    # Constantes de entrada:
    n_prop = 0.8
    gamma = 1.4
    t_t4 = 2460   #[R]
    t_0 = 425      #[R]      
    c_p = 0.24   #[Btu/lbm R]
    h_pr = 18400  #[Btu/lbm]

    #conversões:
    t_0 = t_0/1.8                 #[K]
    c_p  = (c_p  * 4.1868)      #[kJ/kg K]
    c_p  = c_p * 1000              #[J/kg K]
    h_pr = h_pr * 2.326           #[kJ/kg]
    h_pr = h_pr * 1000              #[J/kg]
    t_t4 = t_t4/1.8

    # Tabelas:
    t_t = 0.8
    pi_c_ = np.linspace(2 ,40)  #Limites do eixo X
    lista_empuxo = []
    lista_consumo = []
    plt.figure(figsize = ((12, 6)))
    
    pi_c = 2
    otimo = motor_5.Tau_t_otimo()
    lista_t_t.append(otimo)
   
    for x in range(len(lista_t_t)):
        t_t = lista_t_t[x]
        for i in range(len(pi_c_)):
            pi_c = pi_c_[i]
            lista_empuxo.append(motor_5.Empuxo())
            lista_consumo.append(motor_5.Consumo_especifico())
        if t_t == otimo:
            t_t = format(t_t, ".4f")
            plt.plot(lista_empuxo, lista_consumo, label=r'$\tau_t* = ${}'.format(t_t))
        else:
            plt.plot(lista_empuxo, lista_consumo, label=r'$\tau_t = ${}'.format(t_t))
        lista_empuxo.clear()
        lista_consumo.clear()

    plt.legend(bbox_to_anchor=(1, 1),loc='upper left', borderaxespad=0.)
    plt.xlabel(r'$F/ \dot mt$ ' , fontsize=15)
    plt.ylabel(r'$S$ ' , fontsize=15)
    plt.minorticks_on() # aparece a divisão
    plt.grid()

    # Plot Grafico ----------------------------------------------------------
    plt.show()
    
    print(motor_5.Tau_t_otimo())

    for i in range(len(lista_t_t)):
        t_t = lista_t_t[i]
        print(" ")
        t_t = round(t_t, 3)
        print(f"Razão de temperatura na turbina (tau_t)__ {t_t }")
        print(" ")
        for j in range(len(lista_pi_c)):
            pi_c = lista_pi_c[j]
            print(f"Razão de pressão no combustor (PI_c)____ {pi_c} ")

            dado = format(motor_5.Empuxo(), ".3f")
            print(f"Empuxo (F)______________________________ {dado} [N/(kg s)]")

            dado = format(motor_5.Razao_combustivel_ar(), ".3f")
            print(f"Razão de combustível e ar (f)___________ {dado} ")

            dado = format(motor_5.Consumo_especifico()*1000000, ".3f")
            print(f"Consumo especifico de comustível (S)____ {dado} [(kg s)/(N)]")

            dado = format(motor_5.Eficiencia_termica(), ".3f")
            print(f"Eficiencia termica (n_t)________________ {dado} ")

            dado = format(motor_5.Eficiencia_propulsiva(), ".3f")
            print(f"Eficiencia propulsiva (n_p)_____________ {dado} ")

            dado = format(motor_5.Eficiencia_termica(), ".3f")
            print(f"Eficiencia total (n_o)__________________ {dado} ")

            dado = format(motor_5.C_c(), ".3f")
            print(f"Coeficiente de trabalho (C_c)___________ {dado} ")

            dado = format(motor_5.C_de_prop(), ".3f")
            print(f"Coeficiente de trabalho (C_prop)________ {dado} ")

            dado = format(motor_5.C_de_total(), ".3f")
            print(f"Coeficiente de trabalho total (C_total)_ {dado} ")

            dado = format(otimo, ".4f")
            print(f"Razão de temp total otimo (tau_t*)______ {dado} ")

            print(" ")
        print("**************************************************")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------


