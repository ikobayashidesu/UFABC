import matplotlib.pyplot as plt
import numpy as np

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
        n_t_ = 4
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
    
m_0 = 0.9
pi_f = 4
pi_c = 25

t_0 = 216.7
gamma = 1.4
c_p = 1004.832
h_pr = 42798400
t_t4 = 1666.7

motor = Turbofan_mix_pi_f()
print(motor.Alpha())
print(motor.Razao_combustivel_ar())
print(motor.Razao_de_combustivel_ar_o())
print(motor.Empuxo())
print(motor.Consumo_especifico())

    
#***********************************************************************************************************************************************************