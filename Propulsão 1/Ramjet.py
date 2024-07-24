import matplotlib.pyplot as plt
import numpy as np

#Definindo as constantes de entrada:
# t_0 = 216.7  # [K]
# lambida_ = 1.4
# c_p = 1004   # [J/kg K]
# h_pr = 42800  # [J/kg]

class Ramjet:
    def __init__(self, m_0, t_t4, c_p, lambida_, t_0, h_pr):
        self.m_0 = m_0
        self.t_t4 = t_t4
        self.c_p = c_p
        self.lambida_ = lambida_
        self.t_0 = t_0
        self.h_pr = h_pr

    def Tau_de_r(self):
        tau_r_ = 1 + (((self.lambida_ - 1)/2)*(self.m_0^2))
        self.tau_r = tau_r_
        return tau_r_

    def Tau_de_lambida(self):
        tau_lambida_ = self.t_t4/self.t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_

    def empuxo(self):
        razao_velocidade = self.Raz_vel()
        emp = self.a_0*(razao_velocidade - self.m_0)
        return emp
   
    def Razao_combustivel_ar(self):
        tau_r_ = self.Tau_de_r()
        tau_lambda_ = self.Tau_de_lambida()
        f = ((self.c_p*self.t_0) /self.h_pr)*(tau_lambda_ - tau_r_) 
        self.f_ = f
        return f
    
    def Consumo_especifico(self):
        S = self.Razao_combustivel_ar() / self.empuxo()
        return S
    
motor_1 = Ramjet(2,1600,1004,1.4,216.7,42800)


