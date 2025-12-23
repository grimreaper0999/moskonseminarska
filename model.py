import numpy as np 

def solve_model(T,state):
    mycell_D, mycell_CLK, mycell_MNAND1, mycell_MNAND2, mycell_MNAND3, mycell_MNAND4, mycell_SNAND1, mycell_SNAND2, mycell_Q, mycell_QBAR = state
    dmycell_D = -mycell_D*0
    dmycell_CLK = -mycell_CLK*0
    dmycell_MNAND1 = -mycell_MNAND1*0.1+10*(1)/(1+((mycell_D/5)**3)+((mycell_CLK/5)**3)+((mycell_D/5)**3)*((mycell_CLK/5)**3))
    dmycell_MNAND2 = -mycell_MNAND2*0.1+10*(((mycell_D/5)**3))/(1+((mycell_D/5)**3)+((mycell_CLK/5)**3)+((mycell_D/5)**3)*((mycell_CLK/5)**3))
    dmycell_MNAND3 = -mycell_MNAND3*0.1+10*(1)/(1+((mycell_MNAND1/5)**3)+((mycell_MNAND4/5)**3)+((mycell_MNAND1/5)**3)*((mycell_MNAND4/5)**3))
    dmycell_MNAND4 = -mycell_MNAND4*0.1+10*(1)/(1+((mycell_MNAND2/5)**3)+((mycell_MNAND3/5)**3)+((mycell_MNAND2/5)**3)*((mycell_MNAND3/5)**3))
    dmycell_SNAND1 = -mycell_SNAND1*0.1+10*(((mycell_CLK/5)**3))/(1+((mycell_MNAND3/5)**3)+((mycell_CLK/5)**3)+((mycell_MNAND3/5)**3)*((mycell_CLK/5)**3))
    dmycell_SNAND2 = -mycell_SNAND2*0.1+10*(((mycell_CLK/5)**3))/(1+((mycell_MNAND4/5)**3)+((mycell_CLK/5)**3)+((mycell_MNAND4/5)**3)*((mycell_CLK/5)**3))
    dmycell_Q = -mycell_Q*0.1+10*(1)/(1+((mycell_SNAND1/5)**3)+((mycell_QBAR/5)**3)+((mycell_SNAND1/5)**3)*((mycell_QBAR/5)**3))
    dmycell_QBAR = -mycell_QBAR*0.1+10*(1)/(1+((mycell_SNAND2/5)**3)+((mycell_Q/5)**3)+((mycell_SNAND2/5)**3)*((mycell_Q/5)**3))
    return np.array([dmycell_D, dmycell_CLK, dmycell_MNAND1, dmycell_MNAND2, dmycell_MNAND3, dmycell_MNAND4, dmycell_SNAND1, dmycell_SNAND2, dmycell_Q, dmycell_QBAR])

def solve_model_steady(state):
    return solve_model(0, state)
