import matplotlib.pyplot as plt
import numpy as np

#Definindo as constantes de entrada:
temp_0 = 216.7  # [K]
gamma = 1.4
calor_especifico = 1004   # [J/kg K]
entalpia_total = 42800  # [J/kg]

#Definindo as variaveis de entrada
temp_4 = 1600 # [K]

class Ramjet:
    def __init__(self, t_t4, c_p, lambida_, t_0, h_pr):
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
        raz_2 = mach*raz_1
        self.raz_vel = raz_2
        return raz_2

    def Tau_de_r(self):
        tau_r_ = 1 + (((self.lambida_ - 1)/2)*(mach**2))
        self.tau_r = tau_r_
        return tau_r_

    def Tau_de_lambida(self):
        tau_lambida_ = self.t_t4/self.t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_

    def Empuxo(self):
        razao_velocidade = self.Raz_vel()
        emp = self.A_0()*(razao_velocidade - mach)
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
    
    def Eficiencia_termica(self):
        N_t_ = 1 - (1/self.Tau_de_r())
        self.N_t = N_t_
        return N_t_
    
    def Eficiencia_propulsiva(self):
        N_p_ = 2/(((self.Tau_de_lambida()/self.Tau_de_r())**0.5)+1)
        self.N_p = N_p_
        return N_p_
    
    def Eficiencia_total(self):
        N_0_ = (2*(self.Tau_de_r()-1))/(((self.Tau_de_lambida()*self.Tau_de_r())**0.5)+self.Tau_de_r())
        self.N_0 = N_0_ 
        return N_0_

motor_1 = Ramjet(temp_4,calor_especifico,gamma,temp_0,entalpia_total)
motor_2 = Ramjet(1900,calor_especifico,gamma,temp_0,entalpia_total)
motor_3 = Ramjet(2200,calor_especifico,gamma,temp_0,entalpia_total)

plt.figure(figsize = ((12, 6)))
mach = np.linspace(0 ,7) # grafico de linha (o quanto vai variar o eixo x)

# Grafico 1 ----------------------------------------------------------
plt.subplot(2, 2, 1)
plt.plot(mach, motor_1.Consumo_especifico()*1000, label='1600 K')
plt.plot(mach, motor_2.Consumo_especifico()*1000, label='1900 K')
plt.plot(mach, motor_3.Consumo_especifico()*1000, label='2200 K')
plt.title("Gráfico 01", fontsize = 16)
plt.ylim(40, 100) #limites do eixo Y
plt.legend()
plt.xlabel('M_0')
plt.ylabel('S')
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 2 ----------------------------------------------------------
plt.subplot(2, 2, 2)
plt.plot(mach, motor_1.Razao_combustivel_ar()/1000, label='1600 K')
plt.plot(mach, motor_2.Razao_combustivel_ar()/1000, label='1900 K')
plt.plot(mach, motor_3.Razao_combustivel_ar()/1000, label='2200 K')
plt.title("Gráfico 02", fontsize = 16)
plt.ylim(0, 0.05) #limites do eixo Y
plt.legend()
plt.xlabel('M_0')
plt.ylabel('f')
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 3 ----------------------------------------------------------
plt.subplot(2, 2, 3)
plt.plot(mach, motor_1.Empuxo(), label='1600 K')
plt.plot(mach, motor_2.Empuxo(), label='1900 K')
plt.plot(mach, motor_3.Empuxo(), label='2200 K')
plt.title("Gráfico 03", fontsize = 16)
plt.ylim(0, 1000) #limites do eixo Y
plt.legend()
plt.xlabel('M_0')
plt.ylabel('F/m')
plt.minorticks_on() # aparece a divisão
plt.grid()

# Grafico 4 ----------------------------------------------------------
plt.subplot(2, 2, 4)
plt.plot(mach, motor_1.Eficiencia_propulsiva()*100, label='1600 K')
plt.plot(mach, motor_2.Eficiencia_propulsiva()*100, label='1900 K')
plt.plot(mach, motor_3.Eficiencia_propulsiva()*100, label='2200 K')
plt.plot(mach, motor_1.Eficiencia_total()*100, label='1600 K')
plt.plot(mach, motor_2.Eficiencia_total()*100, label='1900 K')
plt.plot(mach, motor_3.Eficiencia_total()*100, label='2200 K')
plt.title("Gráfico 04", fontsize = 16)
plt.ylim(0, 100) #limites do eixo Y
plt.legend()
plt.xlabel('n')
plt.ylabel('F/m')
plt.minorticks_on() # aparece a divisão
plt.grid()

plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.35)
plt.show()