import matplotlib.pyplot as plt
import numpy as np

#Definindo as constantes de entrada:
temp_0 = 216.7  # [K]
gamma = 1.4
calor_especifico = 1004   # [J/kg K]
entalpia_total = 42800  # [J/kg]

#Definindo as variaveis de entrada
mach = 2
temp_4 = 1600 # [K]

class Ramjet:
    def __init__(self, m_0, t_t4, c_p, lambida_, t_0, h_pr):
        self.m_0 = m_0
        self.t_t4 = t_t4
        self.c_p = c_p
        self.lambida_ = lambida_
        self.t_0 = t_0
        self.h_pr = h_pr

    def R_(self):
        rr = (self.lambida_ - 1)/self.lambida_
        rr = rr * self.c_p
        self.r = rr
        return rr

    def A_0(self):
        rr = self.R_()
        a_0_ = (self.lambida_*rr*self.t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_

    def Raz_vel(self):
        ra = self.Tau_de_r()
        raz = self.Tau_de_lambida()
        raz_1 = (raz/ra)**(0.5)
        raz_2 = self.m_0*raz_1
        self.raz_vel = raz_2
        return raz_2

    def Tau_de_r(self):
        tau_r_ = 1 + (((self.lambida_ - 1)/2)*(self.m_0**2))
        self.tau_r = tau_r_
        return tau_r_

    def Tau_de_lambida(self):
        tau_lambida_ = self.t_t4/self.t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_

    def Empuxo(self):
        razao_velocidade = self.Raz_vel()
        emp = self.A_0()*(razao_velocidade - self.m_0)
        return emp
   
    def Razao_combustivel_ar(self):
        tau_r_ = self.Tau_de_r()
        tau_lambda_ = self.Tau_de_lambida()
        f = ((self.c_p*self.t_0) /self.h_pr)*(tau_lambda_ - tau_r_) 
        self.f_ = f
        return f
    
    def Consumo_especifico(self):
        S = self.Razao_combustivel_ar() / self.Empuxo()
        return S
    
motor_1 = Ramjet(mach,temp_4,calor_especifico,gamma,temp_0,entalpia_total)

print("R:",motor_1.R_())
print("a_0:",motor_1.A_0())
print("t_r:",motor_1.Tau_de_r())
print("t_y:",motor_1.Tau_de_lambida())
print("V_9/a_0:",motor_1.Raz_vel())
print("F/m:",motor_1.Empuxo()) # [N/(kg/s)]
print("f:",motor_1.Razao_combustivel_ar()/1000) 
print("S:",motor_1.Consumo_especifico()*1000) # [(mg/s)/N]
