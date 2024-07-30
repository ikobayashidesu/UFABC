import matplotlib.pyplot as plt
import numpy as np
import time

class Turbojet:
    def __init__(self, tp):
        bora = tp

    def R_(self):
        R_ = ((gamma_-1)/(gamma_))*c_p
        self.R = R_
        return R_

    def A_0(self):
        a_0_ = (gamma_*self.R_()*t_0)**(0.5)
        self.a_0 = a_0_
        return a_0_
    
    def Tau_de_r(self):
        tau_r_ = 1 + (((gamma_ - 1)/(2))*(m_0**2))
        self.tau_r = tau_r_
        return tau_r_
    
    def Tau_de_lambida(self):
        tau_lambida_ = t_t4/t_0
        self.tau_lambida = tau_lambida_
        return tau_lambida_
    
    def Tau_de_c(self):
        tau_c_ = (pi_c)**((gamma_-1)/(gamma_))
        self.tau_c = tau_c_
        return tau_c_
    
    def Tau_de_t(self):
        tau_t_ = 1 - (((self.Tau_de_r())/(self.Tau_de_lambida()))*(self.Tau_de_c()-1))
        self.tau_t = tau_t_
        return tau_t_

    def Razao_de_velocidade(self):
        raz_vel_ = (((2)/(gamma_ -1))*((self.Tau_de_lambida())/(self.Tau_de_r()*self.Tau_de_c()))*((self.Tau_de_c()*self.Tau_de_t()*self.Tau_de_r())-1))**(1/2)
        self.raz_vel = raz_vel_ 
        return raz_vel_ 

    def Empuxo(self):
        empuxo_ = self.A_0()*(self.Razao_de_velocidade()-m_0)
        self.empuxo = empuxo_
        return empuxo_
   
    def Razao_combustivel_ar(self):
        f_ = ((c_p*t_0)/(h_pr))*(self.Tau_de_lambida()-(self.Tau_de_r()*self.Tau_de_c()))
        self.f = f_
        return f_ 
    
    def Consumo_especifico(self):
        s_ = self.Razao_combustivel_ar()/self.Empuxo()
        self.s = s_
        return s_
    
    def Eficiencia_termica(self):
        n_t_ = 1 - ((1)/(self.Tau_de_r()*self.Tau_de_c())) 
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
    
#***********************************************************************************************************************************************************
    
############################################################################################################################################################

print("-------------------------")
print("Digite o numero do motor:")
print("1 - Turbojet")
print("2 - Turbo")
print(" ")
tipo_motor = int(input("-> "))
print("-------------------------")
print(" ")

match tipo_motor:
    case 1:
        print("--Turbojet Selecionado--")
        print(" ")
        print("-------------------------")
        time.sleep(2)
        print(" ")
        print("Definindo as constantes de entrada:")
        print("Temperatura do gás de entrada [K]")
        t_0 = float(input("T_0 = "))
        print(" ")
        print("Valor de gamma:")
        gamma_ = float(input("gamma = "))
        print(" ")
        print("Constante do gás [J/kg K]:")
        c_p = float(input("c_p = "))
        print(" ")
        print("Entalpia total [kJ/kg]:")
        h_pr = float(input("h_pr = "))
        print(" ")
        print("Temperatura de combustão [K]:")
        t_t4 = float(input("T_t4 = "))
        print(" ")
        print(" ")
        print("Definindo as variaveis de entrada:")
        print("Numero de pontos de Mach inicial:")
        numero_mach = int(input("=> "))
        lista_m_0 = []
        lista_m_0_2 = []
        for i in range(numero_mach):
            print(f"faltam {numero_mach - i}:")
            valor = float(input("M_0 = "))
            if valor == 0:
                valor = 0.00000001
            lista_m_0.append(valor)
            lista_m_0_2.append(valor)
        m_0_min = min(lista_m_0)
        m_0_max = max(lista_m_0)
        print(" ")
        print("Razoes PI_c inicial:")
        lista_pi_C = []
        for i in range(numero_mach):
            print(f"faltam {numero_mach - i}:")
            valor_pi = float(input("PI_c = "))
            if valor_pi == 0:
                valor_pi = 0.00000001
            lista_pi_C.append(valor_pi)
        pi_min = min(lista_pi_C)
        pi_max = max(lista_pi_C)
        print(" ")
        print("-------------------------")
        print(" ")


        motor_1 = Turbojet(4)
        plt.figure(figsize = ((12, 6)))
        m_0 = np.linspace(m_0_min ,m_0_max) # grafico de linha (o quanto vai variar o eixo x)

        # Grafico 1 ----------------------------------------------------------
        plt.subplot(2, 2, 1)
        
        for x in range(numero_mach):
            pi_c = lista_pi_C[x]
            plt.plot(m_0, motor_1.Empuxo(), label=r'$\pi_c = ${}'.format(pi_c))
            
        
        plt.ylim(0, 1200) #limites do eixo Y
        plt.legend(bbox_to_anchor=(1.1, 0.2),loc='upper left', borderaxespad=0.)
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$F/ \dot m$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

         # Grafico 2 ----------------------------------------------------------
        plt.subplot(2, 2, 2)
        
        for xd in range(numero_mach):
            pi_c = lista_pi_C[xd]
            plt.plot(m_0, motor_1.Razao_combustivel_ar()/1000)
            
        
        plt.ylim(0, 0.035) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$f$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 3 ----------------------------------------------------------
        plt.subplot(2, 2, 3)

        for xa in range(numero_mach):
            pi_c = lista_pi_C[xa]
            plt.plot(m_0, motor_1.Consumo_especifico()*1000)
            
        
        plt.ylim(0, 100) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$S$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 4 ----------------------------------------------------------
        plt.subplot(2, 2, 4)

        for xb in range(numero_mach):
            pi_c = lista_pi_C[xb]
            plt.plot(m_0, motor_1.Eficiencia_propulsiva(),linestyle=':')

        for xc in range(numero_mach):
            pi_c = lista_pi_C[xc]
            plt.plot(m_0, motor_1.Eficiencia_termica(),linestyle='--')

        for xd in range(numero_mach):
            pi_c = lista_pi_C[xd]
            plt.plot(m_0, motor_1.Eficiencia_total())
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()


        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
        plt.show()

        plt.figure(figsize = ((12, 6)))
        m_0 = np.linspace(m_0_min ,m_0_max) # grafico de linha (o quanto vai variar o eixo x)

        # Grafico EFICIENCIA ----------------------------------------------------------
        # Grafico 1 ----------------------------------------------------------
        plt.subplot(2, 2, 1)

        for xe in range(numero_mach):
            pi_c = lista_pi_C[xe]
            plt.plot(m_0, motor_1.Eficiencia_total(), label=r'$\pi_c = ${}'.format(pi_c))
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.legend(bbox_to_anchor=(1.1, 0.2),loc='upper left', borderaxespad=0.)
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta_t$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

         # Grafico 2 ----------------------------------------------------------
        plt.subplot(2, 2, 2)

        for xf in range(numero_mach):
            pi_c = lista_pi_C[xf]
            plt.plot(m_0, motor_1.Eficiencia_propulsiva(),linestyle=':')
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta_p$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 3 ----------------------------------------------------------
        plt.subplot(2, 2, 3)
        for xg in range(numero_mach):
            pi_c = lista_pi_C[xg]
            plt.plot(m_0, motor_1.Eficiencia_termica(),linestyle='--')
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta_o$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 4 ----------------------------------------------------------
        plt.subplot(2, 2, 4)
        for xh in range(numero_mach):
            pi_c = lista_pi_C[xh]
            plt.plot(m_0, motor_1.Eficiencia_propulsiva(),linestyle=':')

        for xi in range(numero_mach):
            pi_c = lista_pi_C[xi]
            plt.plot(m_0, motor_1.Eficiencia_termica(),linestyle='--')

        for xj in range(numero_mach):
            pi_c = lista_pi_C[xj]
            plt.plot(m_0, motor_1.Eficiencia_total())
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()


        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
        plt.show()
        
        mx = 0
        for mx in range(numero_mach):
            m_0 = lista_m_0_2[mx]
            print(f"Número de Mach = {m_0}")
            for xpi in range(numero_mach):
                pi_c = lista_pi_C[xpi]
                print(f"Razão de pressão do compressor = {pi_c}")
                resultado = format(motor_1.Empuxo(), ".3f")
                print(f"Empuxo (F) = {resultado} N/(kg/s)")

                resultado = format(motor_1.Razao_combustivel_ar()/1000, ".3f")
                print(f"Razão combustivel/ar (f) = {resultado}")

                resultado = format(motor_1.Consumo_especifico()*1000, ".3f")
                print(f"Consumo especifico de combustivel (S) = {resultado} (mg/2)/N ")

                resultado = format(motor_1.Eficiencia_termica(), ".3f")
                print(f"Eficiencia termica (n_t) = {resultado}")
                resultado = format(motor_1.Eficiencia_propulsiva(), ".3f")
                print(f"Eficiencia propulsiva (n_p) = {resultado}")
                resultado = format(motor_1.Eficiencia_total(), ".3f")
                print(f"Eficiencia total (n_o) = {resultado}")
                
                print(" ")
            mx = mx + 1
            print("-----------------------------------------")
    
    case 2:
        print("--Turbojet Selecionado--")
        print(" ")
        print("-------------------------")
        time.sleep(2)
        print(" ")
        print("Definindo as constantes de entrada:")
        print("Temperatura do gás de entrada [K]")
        t_0 = float(input("T_0 = "))
        print(" ")
        print("Valor de gamma:")
        gamma_ = float(input("gamma = "))
        print(" ")
        print("Constante do gás [J/kg K]:")
        c_p = float(input("c_p = "))
        print(" ")
        print("Entalpia total [kJ/kg]:")
        h_pr = float(input("h_pr = "))
        print(" ")
        print("Temperatura de combustão [K]:")
        t_t4 = float(input("T_t4 = "))
        print(" ")
        print(" ")
        print("Definindo as variaveis de entrada:")
        print("Numero de pontos de Mach inicial:")
        numero_mach = int(input("=> "))
        lista_m_0 = []
        lista_m_0_2 = []
        for i in range(numero_mach):
            print(f"faltam {numero_mach - i}:")
            valor = float(input("M_0 = "))
            if valor == 0:
                valor = 0.00000001
            lista_m_0.append(valor)
            lista_m_0_2.append(valor)
        m_0_min = min(lista_m_0)
        m_0_max = max(lista_m_0)
        print(" ")
        print("Razoes PI_c inicial:")
        lista_pi_C = []
        for i in range(numero_mach):
            print(f"faltam {numero_mach - i}:")
            valor_pi = float(input("PI_c = "))
            if valor_pi == 0:
                valor_pi = 0.00000001
            lista_pi_C.append(valor_pi)
        pi_min = min(lista_pi_C)
        pi_max = max(lista_pi_C)
        print(" ")
        print("-------------------------")
        print(" ")


        motor_1 = Turbojet(4)
        plt.figure(figsize = ((12, 6)))
        m_0 = np.linspace(m_0_min ,m_0_max) # grafico de linha (o quanto vai variar o eixo x)

        # Grafico 1 ----------------------------------------------------------
        plt.subplot(2, 2, 1)
        
        for x in range(numero_mach):
            pi_c = lista_pi_C[x]
            plt.plot(m_0, motor_1.Empuxo(), label=r'$\pi_c = ${}'.format(pi_c))
            
        
        plt.ylim(0, 1200) #limites do eixo Y
        plt.legend(bbox_to_anchor=(1.1, 0.2),loc='upper left', borderaxespad=0.)
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$F/ \dot m$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

         # Grafico 2 ----------------------------------------------------------
        plt.subplot(2, 2, 2)
        
        for xd in range(numero_mach):
            pi_c = lista_pi_C[xd]
            plt.plot(m_0, motor_1.Razao_combustivel_ar()/1000)
            
        
        plt.ylim(0, 0.035) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$f$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 3 ----------------------------------------------------------
        plt.subplot(2, 2, 3)

        for xa in range(numero_mach):
            pi_c = lista_pi_C[xa]
            plt.plot(m_0, motor_1.Consumo_especifico()*1000)
            
        
        plt.ylim(0, 100) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$S$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()

        # Grafico 4 ----------------------------------------------------------
        plt.subplot(2, 2, 4)

        for xb in range(numero_mach):
            pi_c = lista_pi_C[xb]
            plt.plot(m_0, motor_1.Eficiencia_propulsiva(),linestyle=':')

        for xc in range(numero_mach):
            pi_c = lista_pi_C[xc]
            plt.plot(m_0, motor_1.Eficiencia_termica(),linestyle='--')

        for xd in range(numero_mach):
            pi_c = lista_pi_C[xd]
            plt.plot(m_0, motor_1.Eficiencia_total())
            
        
        plt.ylim(0, 1) #limites do eixo Y
        plt.xlabel(r'$M_0$ ' , fontsize=15)
        plt.ylabel(r'$\eta$ ' , fontsize=15)
        plt.minorticks_on() # aparece a divisão
        plt.grid()


        plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.7, hspace=0.35)
        plt.show()

            

        



 

    case _:
        print("-- Nâo existe essa opção --")

      